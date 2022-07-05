

class TLEException(Exception):
    """Exception raised for errors in parsing a TLE.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class TokenNumberException(TLEException):
    """Exception raised for an invalid number of tokens in a TLE line.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class TokenLengthException(TLEException):
    """Exception raised for an invalid number of characters in a TLE token.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)
