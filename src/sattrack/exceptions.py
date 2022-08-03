

class TLEException(Exception):
    """
    Exception raised for errors in parsing a TLE.

    Attributes
        message: Explanation of the error.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class LineNumberException(TLEException):
    """
    Exception raised for an invalid number of lines in a TLE.

    Attributes:
        message: Explanation of the error.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class ChecksumException(TLEException):
    """
    Exception raised for an invalid checksum for a TLE line.

    Attributes:
        message: Explanation of the error.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class TokenNumberException(TLEException):
    """
    Exception raised for an invalid number of tokens in a TLE line.

    Attributes:
        message: Explanation of the error.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class TokenLengthException(TLEException):
    """
    Exception raised for an invalid number of characters in a TLE token.

    Attributes:
        message: Explanation of the error.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class NoPassException(Exception):
    """
    Exception raised for when an overhead satellite pass does not occur.

    Attributes
        message: Explanation of the error.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)
