from copy import copy
from math import cos, sin

from pyevspace import EVector, EMatrix, transpose

from sattrack.rotation.order import Axis, EulerOrder

__all__ = ('EulerAngles', 'getMatrix', 'getEulerMatrix', 'getMatrixFromTo', 'rotateAxisTo', 'rotateAxisFrom',
           'rotateMatrixTo', 'rotateMatrixFrom', 'rotateOrderTo', 'rotateOrderFrom', 'rotateFromTo', 'ReferenceFrame',
           'rotateToThenOffset', 'undoRotateToThenOffset')


def _xRotation(angle: float) -> EMatrix:
    angleCos = cos(angle)
    angleSin = sin(angle)
    return EMatrix(
            [1.0, 0.0, 0.0],
            [0.0, angleCos, -angleSin],
            [0.0, angleSin, angleCos]
    )


def _yRotation(angle: float) -> EMatrix:
    angleCos = cos(angle)
    angleSin = sin(angle)
    return EMatrix(
            [angleCos, 0.0, angleSin],
            [0.0, 1.0, 0.0],
            [-angleSin, 0.0, angleCos]
    )


def _zRotation(angle: float) -> EMatrix:
    angleCos = cos(angle)
    angleSin = sin(angle)
    return EMatrix(
            [angleCos, -angleSin, 0.0],
            [angleSin, angleCos, 0.0],
            [0.0, 0.0, 1.0]
    )


class EulerAngles:
    """
    A container class to hold the rotation angle values of an Euler rotation.

    The angle values correspond to the axes in the EulerOrder object, i.e. the index of the angle should match the index
    of the axis in the corresponding EulerOrder. All angles are measured in radians.
    """

    def __init__(self, alpha: float, beta: float, gamma: float):
        """
        Initializes the angles to the parameters given.
        Args:
            alpha: The angle of the first rotation in radians.
            beta: The angle of the second rotation in radians.
            gamma: The angle of the third rotation in radians.
        """

        self._angles = [alpha, beta, gamma]

    def __getitem__(self, index: int) -> float:
        """
        Returns the specified rotation angle.

        Args:
            index: The index of the desired angle.

        Returns:
            The rotation angle in radians.

        Raises:
            ValueError: If the index is not in the range (0, 2).
        """
        if index < 0 or index > 2:
            raise ValueError("Index value out of range.")
        return self._angles[index]

    def __setitem__(self, index: int, value: float) -> None:
        """
        Sets the specified rotation angle.

        Args:
            index: The index of the angle to change.
            value: The value to set the angle in radians.

        Raises:
            ValueError: If the index is not in the range (0, 2).
        """
        if index < 0 or index > 2:
            raise ValueError("Index value out of range.")
        self._angles[index] = value

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

    def __reduce__(self):
        return self.__class__, (self._angles[0], self._angles[1], self._angles[2])


# purposefully not documenting these since I intend to implement this in C.
# is an extrinsic rotation
def getMatrix(axis: Axis, angle: float) -> EMatrix:
    """Creates a rotation matrix for rotating extrinsically around a single axis.
    Parameters:
    axis:   The axis to rotate around, an Axis enumerated type.
    angle:  The angle to rotate by, measured in radians."""
    if axis == Axis.X_AXIS:
        return _xRotation(angle)
    elif axis == Axis.Y_AXIS:
        return _yRotation(angle)
    elif axis == Axis.Z_AXIS:
        return _zRotation(angle)
    else:
        return EMatrix.I


# is made up of intrinsic rotations
def getEulerMatrix(order: EulerOrder, angles: EulerAngles) -> EMatrix:
    """Creates a rotation matrix based on the Euler rotation order and the respective rotation angles.
    Parameters:
    order:  An order.EulerOrder object in order.Order corresponding to the specific axes and their order of rotation.
    angles: A rotation.EulerAngles object containing the angles corresponding to the axis rotations."""
    rtn = EMatrix.I
    for axis, angle in zip(order, angles):
        rtn = rtn * getMatrix(axis, angle)
    return rtn


def getMatrixFromTo(orderFrom: EulerOrder, angsFrom: EulerAngles, orderTo: EulerOrder, angsTo: EulerAngles) -> EMatrix:
    """Creates a rotation matrix describing a rotation between two reference frames.
    Parameters:
    orderFrom:  An order.EulerOrder object in order.Order corresponding to the reference frame rotating from.
    angsFrom:   A rotation.EulerAngles object containing the angles corresponding to the axis rotations for the
                reference frame rotating from.
    orderTo:    An order.EulerOrder object in order.Order corresponding to the reference frame rotating to.
    angsTo:     A rotation.EulerAngles object containing the angles corresponding to the axis rotations for the
                reference frame rotating to.
    """
    return transpose(getEulerMatrix(orderTo, angsTo)) * getEulerMatrix(orderFrom, angsFrom)


def rotateAxisTo(axis: Axis, angle: float, vector: EVector) -> EVector:
    """Returns a vector rotated from an inertial reference frame 'to' the reference frame
    represented by the single axis rotation.
    Parameters:
    axis:   One of the order.Axis enumerated types.
    angle:  The angle of rotation measured in radians.
    vector: The vector to rotate."""
    return transpose(getMatrix(axis, angle)) * vector


def rotateAxisFrom(axis: Axis, angle: float, vector: EVector) -> EVector:
    """Returns a vector rotated to an inertial reference frame 'from' the reference frame
    represented by the single axis rotation.
    Parameters:
    axis:   One of the order.Axis enumerated types.
    angle:  The angle of rotation measured in radians.
    vector: The vector to rotate."""
    return getMatrix(axis, angle) * vector


def rotateMatrixTo(matrix: EMatrix, vector: EVector) -> EVector:
    """Returns a vector rotated from an inertial reference frame 'to' the reference frame
    represented by the rotation matrix.
    Parameters:
    matrix: The rotation matrix.
    vector: The vector to rotate."""
    return transpose(matrix) * vector


def rotateMatrixFrom(matrix: EMatrix, vector: EVector) -> EVector:
    """Returns a vector rotated to an inertial reference frame 'from' the reference frame
    represented by the rotation matrix.
    Parameters:
    matrix: The rotation matrix.
    vector: The vector to rotate."""
    return matrix * vector


def rotateOrderTo(order: EulerOrder, angles: EulerAngles, vector: EVector) -> EVector:
    """Returns a vector rotated from an inertial reference frame 'to' the reference frame
    represented by the euler rotation.
    Parameters:
    order:  An order.EulerOrder object found in order.Order representing the rotation order.
    angles: A rotation.EulerAngles object containing the angles corresponding to the axis rotations.
    vector: The vector to rotate."""
    return transpose(getEulerMatrix(order, angles)) * vector


def rotateOrderFrom(order: EulerOrder, angles: EulerAngles, vector: EVector) -> EVector:
    """Returns a vector rotated to an inertial reference frame 'from' the reference frame
    represented by the euler rotation.
    Parameters:
    order:  An order.EulerOrder object found in order.Order representing the rotation order.
    angles: A rotation.EulerAngles object containing the angles corresponding to the axis rotations.
    vector: The vector to rotate."""
    return getEulerMatrix(order, angles) * vector


def rotateFromTo(orderFrom: EulerOrder, angsFrom: EulerAngles, orderTo: EulerOrder, angsTo: EulerAngles,
                 vector: EVector) -> EVector:
    """Creates a rotation matrix corresponding to a rotation from a reference frame to another.
    Parameters:
    orderFrom:  An order.EulerOrder object found in order.Order representing the reference frame rotating from.
    angsFrom:   A rotation.EulerAngles object containing the angles corresponding to the axis rotations """
    return getMatrixFromTo(orderFrom, angsFrom, orderTo, angsTo) * vector


class ReferenceFrame:
    """
    Class that represents a rotation between reference frames.

    Each of the rotation components are stored in and are accessible through the class. The angles can be adjusted, as
    they naturally tend to do for a non-inertial reference frame, however the rotation order is not expected to change,
    so they are not mutable after instantiation. The benefit of the class is the ability to conveniently store all
    objects relating to a rotation, as well as the built-in methods for rotating EVector's between the reference frames.
    """

    def __init__(self, order: EulerOrder, angles: EulerAngles):
        """Initializes the internal matrix based on the order and angle arguments."""
        self._order = order
        self._angles = angles
        self._matrix = getEulerMatrix(order, angles)

    def setAngles(self, angles: EulerAngles):
        """Sets the rotation angles and updates the internal rotation matrix."""
        self._angles = angles
        self._matrix = getEulerMatrix(self._order, angles)

    def getAngles(self) -> EulerAngles:
        """Returns a copy of the EulerAngles for this rotation."""
        return copy(self._angles)

    def getOrder(self) -> EulerOrder:
        """Returns the EulerOrder that describes this rotation."""
        return self._order

    def getMatrix(self) -> EMatrix:
        """Returns the rotation matrix that corresponds to this rotation."""
        return self._matrix.copy()

    def RotateTo(self, vector: EVector) -> EVector:
        """
        Rotates a vector from inertial reference frame coordinates, to the rotated reference frame coordinates.

        For rotating to a non-inertial reference frame, use the RotateToFrame() class method.

        Args:
            vector: Vector to be rotated.

        Returns:
            The vector rotated to the reference frame represented by the instance.
        """

        return transpose(self._matrix) * vector

    def RotateFrom(self, vector: EVector) -> EVector:
        """
        Rotates a vector from the rotated reference frame's coordinates, to inertial reference frame coordinates.

        For rotating from a non-inertial reference frame, use the RotateFromFrame() class method.

        Args:
            vector: Vector to be rotated.

        Returns:
            The vector rotated from the reference frame represented by the instance.
        """

        return self._matrix * vector

    def RotateToFrame(self, refFrame, vector: EVector) -> EVector:
        """
        Rotates a vector from this reference frame's coordinates, to another reference frame's coordinates.

        For rotating to an inertial reference frame, use the RotateTo() class method.

        Args:
            refFrame: Non-inertial ReferenceFrame object.
            vector: Vector to be rotated

        Returns:
            The vector rotated to the ReferenceFrame parameter's reference frame.
        """

        return transpose(refFrame.getMatrix()) * self._matrix * vector

    def RotateFromFrame(self, refFrame, vector: EVector) -> EVector:
        """
        Rotates a vector from a reference frame's coordinates, to this reference frame's coordinates.

        For rotating from an inertial reference frame, use the RotateFrame() class method.

        Args:
            refFrame: Non-inertial ReferenceFrame object.
            vector: Vector to be rotated.

        Returns:
            The vector rotated from the ReferenceFrame parameter's reference frame.
        """

        return transpose(self._matrix) * refFrame.getMatrix() * vector


def rotateToThenOffset(rotation: EMatrix, offset: EVector, original: EVector) -> EVector:
    """
    Rotates a vector to a rotated and offset reference frame.

    Args:
        rotation: Rotation matrix of offset reference frame.
        offset: Offset vector between reference frame origins.
        original: Original vector to be rotated.

    Returns:
        The original vector in the rotated and offset reference frame.
    """

    rotatedOriginal = transpose(rotation) * original
    rotatedOffset = transpose(rotation) * offset
    return rotatedOriginal - rotatedOffset


def undoRotateToThenOffset(rotation: EMatrix, offset: EVector, rotated: EVector) -> EVector:
    """
    Un-offsets a vector and then rotates it back into its original reference frame. This effectively undoes what the
    rotateToThenOffset method does to its 'original' parameter.

    Args:
        rotation: Rotation matrix of offset reference frame.
        offset: Offset vector between reference frame origins.
        rotated: Rotated and offset vector to be returned to original reference frame.

    Returns:
        The original vector.
    """

    rotatedOffset = transpose(rotation) * offset
    rotatedOriginal = rotated + rotatedOffset
    return rotation * rotatedOriginal
