import glob
import configparser
from .errors import hffError

try:
    import astromorph
    _HAS_ASTROMORPH = True
except ImportError:
    _HAS_ASTROMORPH = False


if _HAS_ASTROMORPH is True:
    from astromorph.lensing import LensingModel,stack_models


class HFFData:

    cluster_coords = {"macs0717":(109.43863048556054, 109.33323737580649, 37.708234125318704, 37.791567413408224),\
                      "abell2744":(3.6349786790793948, 3.531611340899466, -30.439247968007468, -30.336748088394312),\
                      "abell370":(40.00302643255294, 39.92299569762585, -1.6432928825364999, -1.5332929431505387),\
                      "macs1149":(177.44326387785338, 177.34591790002193, 22.356309256746993, 22.446309208248984),\
                      "abell1063":(342.23815044338255, 342.12452654293924, -44.57507016811651, -44.48856837273851),\
                      "macs0416":(64.08056990070776, 63.989300170823576, -24.114107728713556, -24.030773153403377)}


    def __init__(self,dataFolder):
        self.folder = dataFolder
        self.lensing_data = {}
        return

    def get_cluster_from_coords(self,coords):
        ra = coords.ra.value
        # dec = coords.dec.value
        for cluster,limits in HFFData.cluster_coords.items():
            if (ra>limits[1]) and (ra<limits[0]):
                return cluster
        raise hffError(f"Coordinates {ra},{dec} outside cluster boundaries.")
        return None

    def find_all_models(self,rejectList=[]):
        allFiles = glob.glob(f"{self.folder}/*/v*/*.fits")
        roots = []
        for f in allFiles:
            reject = False
            fname = f.split("/")[-1]
            rootName = "_".join(fname.split("_")[:6])
            for model in rejectList:
                if model in rootName:
                    reject = True
                    break
            if reject:
                continue
            if not rootName in roots:
                roots.append(rootName)
        return roots

    def load_data(self,cluster,rejectModels=[]):
        lens_models_folder = f"{self.folder}/{cluster}"
        print(lens_models_folder)
        lensModelNames = self.find_all_models(rejectModels)
        print(lensModelNames)
        nModels = len(lensModelNames)
        ConfigFile = configparser.ConfigParser()
        ConfigFile.read("%s/models.cfg"%(lens_models_folder))
        modelsLensing = []
        for i in range(nModels):
            modelName = lensModelNames[i].split("_")[-2]
            version = lensModelNames[i].split("_")[-1]
            redshiftLens = float(ConfigFile[modelName]["redshift"])
            resolution = ConfigFile[modelName]["resolution"]
            if len(resolution.split(","))==1:
                resolutionModel = float(ConfigFile[modelName]["resolution"])
            else:
                resolutionModel = float(ConfigFile[modelName]["resolution"].split(",")[0])
            lensModel = LensingModel(redshiftLens,resolutionModel)
            lensModel.set_lensing_data(filename=f"{lens_models_folder}/{modelName}/{version}/{lensModelNames[i]}")
            modelExtent = lensModel.get_image_box_coordinates()
            modelsLensing.append(lensModel)

        self.lensing_data[cluster] = {}
        self.lensing_data[cluster]["models"] = modelsLensing
        self.lensing_data[cluster]["redshift"] = redshiftLens
        return modelsLensing


    def get_lensing_pars(self,coords,size,redshift,pixScaleLensStack = 0.2,rejectModels=[]):
        if _HAS_ASTROMORPH is False:
            raise hffError("Need to have astromorph package installed. Check https://github.com/AfonsoV/astromorph.")
        cluster = self.get_cluster_from_coords(coords)
        if not cluster in self.lensing_data.keys():
            modelData = self.load_data(cluster,rejectModels=rejectModels)
        else:
            modelData = self.lensing_data[cluster]

        stackRA = (coords.ra.value-size/(2*3600.),coords.ra.value+size/(2*3600.))
        stackDEC = (coords.dec.value-size/(2*3600.),coords.dec.value+size/(2*3600.))

        print("\tStacking Lens Models...")
        modelStack = stack_models(lensModels,stackRA,stackDEC,\
                                  scale=pixScaleLensStack/3600.0,\
                                  modelbbox=(1.1*size,coords))
        ExtentModel = stackRA+stackDEC

        LensParsGalaxyFull = []
        for i in range(3):
            LensingModelStacked = LensingModel(lensRedshift,pixScaleLensStack)
            LensingModelStacked.set_lensing_data(kappa=modelStack[i,0,:,:],gamma=modelStack[i,1,:,:],\
                                                 xdeflect=modelStack[i,2,:,:],ydeflect=modelStack[i,3,:,:],\
                                                 extent=ExtentModel)
            LensingModelStacked.set_model_at_z(redshift)
            LensingModelStacked.compute_shear_angle(redshift)
            LensParsGalaxy = LensingModelStacked.get_lensing_parameters_at_position(coords)
            LensParsGalaxyFull.append(LensParsGalaxy)

        LensParsGalaxyFull = np.asarray(LensParsGalaxyFull)
        errorsLensPars = LensParsGalaxyFull[1:]-LensParsGalaxyFull[:-1]
        LensPars = LensParsGalaxyFull[1,:]
        print(errorsLensPars.shape)
        return LensPars,errorsLensPars

    def get_magnification(self,coords,size,redshift,pixScaleLensStack = 0.2,rejectModels=[]):
        pars,errors = self.get_lensing_pars(coords,size,redshift,\
                              pixScaleLensStack = pixScaleLensStack,\
                              rejectModels=rejectModels)

        mu = pars[2]
        mu_err = errors[2,:]
        return mu,mu_err
