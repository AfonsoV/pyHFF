from setuptools import setup

setup(
  name = "pyHFF",
  packages = ["hff"],
  version = "0.1.0",
  description = "A python set of utilities to deal with HFF data",
  author = "Bruno Ribeiro",
  author_email = "brunorlr@gmail.com",
  url = "https://github.com/AfonsoV/pyHFF", # use the URL to the github repo
  # download_url = "https://github.com/afonsov/astromorph/archive/0.1.tar.gz",
  keywords = ["astronomy", "galaxies", "HFF"],
  classifiers = ["Development Status :: 2 - Pre-Alpha",\
                 "Programming Language :: Python :: 3",\
                 "Topic :: Scientific/Engineering :: Astronomy",\
                 "Intended Audience :: Science/Research",\
                 "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"],
  license="GPLv3",
  include_package_data=False,
  install_requires=[]
)
