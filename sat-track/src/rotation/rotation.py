from order import Axis, Order, EulerOrder
from pyevspace import EVector, EMatrix, transpose
from math import cos, sin, radians
from copy import deepcopy

# todo: change this when we can create matrices from lists
def _xrotation(angle: float) -> EMatrix:
    rtn = EMatrix()
    rtn.set(0, 0, 1.0)
    anglerad = radians(angle)
    anglecos = cos(anglerad)
    anglesin = sin(anglerad)
    rtn.set(1, 1, anglecos)
    rtn.set(2, 2, anglecos)
    rtn.set(1, 2, -anglesin)
    rtn.set(2, 1, anglesin)
    return rtn

def _yrotation(angle: float) -> EMatrix:
    rtn = EMatrix()
    rtn.set(1, 1, 1.0)
    anglerad = radians(angle)
    anglecos = cos(anglerad)
    anglesin = sin(anglerad)
    rtn.set(0, 0, anglecos)
    rtn.set(2, 2, anglecos)
    rtn.set(2, 0, -anglesin)
    rtn.set(0, 2, anglesin)
    return rtn

def _zrotation(angle: float) -> EMatrix:
    rtn = EMatrix()
    rtn.set(2, 2, 1.0)
    anglerad = radians(angle)
    anglecos = cos(anglerad)
    anglesin = sin(anglerad)
    rtn.set(0, 0, anglecos)
    rtn.set(1, 1, anglecos)
    rtn.set(0, 1, -anglesin)
    rtn.set(1, 0, anglesin)
    return rtn

class EulerAngles:

    def __init__(self, alpha: float, beta: float, gamma: float):
        self._angles = [alpha, beta, gamma]

    def __getitem__(self, i: int) -> float:
        if (i < 0 or i > 2):
            raise ValueError("Index value out of range.")
        return self._angles[i]

    def __setitem__(self, i:int, val: float) -> None:
        if (i < 0 or i > 2):
            raise ValueError("Index value out of range.")
        self._angles[i] = val

    def __str__(self) -> str:
        return str([i for i in self._angles])

    def __iter__(self):
        self._n = 0
        return self

    def __next__(self):
        if self._n <= 2:
            temp = self._angles[self._n]
            self._n += 1
            return temp
        else:
            raise StopIteration

# is an extrinsic rotation
def getMatrix(axis: Axis, angle: float) -> EMatrix:
    '''Creates a rotation matrix for rotating extrinsically around a single axis.
    Parameters:
    axis:   The axis to rotate around, an Axis enumerated type.
    angle:  The angle to rotate by, measured in degrees.'''
    if axis == Axis.X_AXIS:
        return _xrotation(angle)
    elif axis == Axis.Y_AXIS:
        return _yrotation(angle)
    elif axis == Axis.Z_AXIS:
        return _zrotation(angle)
    else:
        return EMatrix.I

# is made up of intrinsic rotations
def getEulerMatrix(order: EulerOrder, angles: EulerAngles) -> EMatrix:
    '''Creates a rotation matrix based on the Euler rotation order and the respective rotation angles.
    Parameters:
    order:  An order.EulerOrder object in order.Order cooresponding to the specific axes and their order of rotation.
    angles: A rotation.EulerAngles object containing the angles cooresponding to the axis rotations.'''
    rtn = EMatrix.I
    for axis, angle in zip(order, angles):
        rtn = rtn @ getMatrix(axis, angle)
    return rtn;

def getMatrixFromTo(orderFrom: EulerOrder, angsFrom: EulerAngles, orderTo: EulerOrder, angsTo: EulerAngles) -> EVector:
    '''Creates a rotation matrix describing a rotation between two reference frames.
    Parameters:
    orderFrom:  An order.EulerOrder object in order.Order cooresponding to the reference frame rotating from.
    angsFrom:   A rotation.EulerAngles object containing the angles cooresponding to the axis rotations for the
                reference frame rotating from.
    orderTo:    An order.EulerOrder object in order.Order cooresponding to the reference frame rotating to.
    angsTo:     A rotation.EulerAngles object containing the angles cooresponding to the axis rotatinos for the
                reference frame rotating to.
    '''
    return transpose(getEulerMatrix(orderTo, angsTo)) @ getEulerMatrix(orderFrom, angsFrom)

def rotateAxisTo(axis: Axis, angle: float, vector: EVector) -> EVector:
    '''Returns a vector rotated from an inertial reference frame 'to' the reference frame
    represented by the single axis rotation.
    Parameters:
    axis:   One of the order.Axis enumerated types.
    angle:  The angle of rotation measured in degrees.
    vector: The vector to rotate.'''
    return transpose(getMatrix(axis, angle)) @ vector

def rotateAxisFrom(axis: Axis, angle: float, vector: EVector) -> EVector:
    '''Returns a vector rotated to an inertial reference frame 'from' the reference frame
    represented by the single axis rotation.
    Parameters:
    axis:   One of the order.Axis enumerated types.
    angle:  The angle of rotation measured in degrees.
    vector: The vector to rotate.'''
    return getMatrix(axis, angle) @ vector

def rotateMatrixTo(matrix: EMatrix, vector: EVector) -> EVector:
    '''Returns a vector rotated from an inertial reference frame 'to' the reference frame
    represented by the rotation matrix.
    Parameters:
    matrix: The rotation matrix.
    vector: The vector to rotate.'''
    return transpose(matrix) @ vector

def rotateMatrixFrom(matrix: EMatrix, vector: EVector) -> EVector:
    '''Returns a vector rotated to an inertial reference frame 'from' the reference frame
    represented by the rotation matrix.
    Parameters:
    matrix: The rotation matrix.
    vector: The vector to rotate.'''
    return matrix @ vector

def rotateOrderTo(order: EulerOrder, angles: EulerAngles, vector: EVector) -> EVector:
    '''Returns a vector rotated from an inertial reference frame 'to' the reference frame
    represented by the euler rotation.
    Parameters:
    order:  An order.EulerOrder object found in order.Order representing the rotation order.
    angles: A rotation.EulerAngles object containing the angles cooresponding to the axis rotations.
    vector: The vector to rotate.'''
    return transpose(getEulerMatrix(order, angles)) @ vector

def rotateOrderFrom(order: EulerOrder, angles: EulerAngles, vector: EVector) -> EVector:
    '''Returns a vector rotated to an inertial reference frame 'from' the reference frame 
    represented by the euler rotation.
    Parameters:
    order:  An order.EulerOrder object found in order.Order representing the rotation order.
    angles: A rotation.EulerAngles object containing the angles coordsponding to the axis rotations.
    vector: The vector to rotate.'''
    return getEulerMatrix(order, angles) @ vector

def rotateFromTo(orderFrom: EulerOrder, angsFrom: EulerAngles, orderTo: EulerOrder, angsTo: EulerAngles, vector: EVector) -> EVector:
    '''Creates a rotation matrix cooresponding to a rotation from a reference frame to another.
    Parameters:
    orderFrom:  An order.EulerOrder object found in order.Order representing the reference frame rotating from.
    angsFrom:   A rotation.EulerAngles object containing the angles cooresponding to the axis rotations '''
    return getMatrixFromTo(orderFrom, angsFrom, orderTo, angsTo) @ vector

class Rotation:
    '''Rotation object that contains a rotation order, and the angles of the rotations.
    This object also contains an internal rotation matrix representing the rotation.
    The angles can be adjusted, as they naturally tend to do for a non-inertial reference frame,
    however the rotation order is not expected to change. Therefore, if you wish to represent a 
    rotation that has a different order, you should create a new Rotation object altogether.'''
    
    def __init__(self, order: EulerOrder, angles: EulerAngles):
        self._order = order
        self._angles = angles
        self._matrix = getEulerMatrix(order, angles)

    def setAngles(self, angles: EulerAngles):
        self._angles = angles
        self._matrix = getEulerMatrix(self._order, angles)

    def getAngles(self) -> EulerAngles:
        return deepcopy(self._angles)

    def getOrder(self) -> EulerOrder:
        return self._order

    def getMatrix(self) -> EMatrix:
        return deepcopy(self._matrix)

    # rotate vector to the reference frame represented by the rotation
    def RotateTo(self, vector: EVector) -> EVector:
        return transpose(self._matrix) @ vector

    # rotate a vector from the reference frame represented by the rotation
    def RotateFrom(self, vector: EVector) -> EVector:
        return self._matrix @ vector

    def RotateToFrame(self, rotation, vector: EVector) -> EVector:
        return transpose(rotation._matrix) @ self._matrix @ vector

    def RotateFromFrame(self, rotation, vector: EVector) -> EVector:
        return transpose(self._matrix) @ rotation._matrix @ vector


def rotateToThenOffset(rotation: EMatrix, offset: EVector, original: EVector) -> EVector:
    rotatedOriginal = transpose(rotation) @ original
    rotatedOffset = transpose(rotation) @ offset
    return rotatedOriginal - rotatedOffset
    