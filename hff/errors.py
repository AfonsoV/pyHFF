

class hffException(Exception):
    """Base class for exceptions in this module."""
    pass

class hffError(hffException):
    """Exception raised for general errors.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
