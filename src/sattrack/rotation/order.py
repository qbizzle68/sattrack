from enum import Enum


__all__ = ('Axis', 'EulerOrder', 'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX', 'XYX', 'XZX', 'YXY', 'YZY', 'ZXZ', 'ZYZ')


# enum for general axis directions, empty_axis has been useful before (2 rotations)
class Axis(Enum):
    """
    Enumeration to differentiate the three axis directions.

    The enumeration includes an empty axis, which can be used when an Axis type is needed programmatically, but
    logically the type is not significant. For example an Euler rotation order of only two axes instead of three can
    be described with the last rotation Axis type being EMPTY_AXIS.
    """

    X_AXIS = 0
    Y_AXIS = 1
    Z_AXIS = 2
    EMPTY_AXIS = 3


class EulerOrder:
    """
    A container class used describe the order of an Euler rotation.

    The object is instantiated with three axes, representing the three axes, in order, of the Euler rotation. The
    __iter__ and __next__ methods have been implemented so an instance can be iterated over, and the __getitem__ method
    can be used to retrieve any desired axis.
    """

    def __init__(self, axis1: Axis, axis2: Axis, axis3: Axis):
        """
        Initializes an instance whose axis arguments are the axes of an Euler rotation.

        Args:
            axis1: The first axis of the Euler rotation.
            axis2: The second axis of the Euler rotation.
            axis3: The third axis of the Euler rotation.
        """
        self._rotation_order = (axis1, axis2, axis3)

    def __getitem__(self, index: int) -> Axis:
        """
        Returns the specified axis of the rotation.

        Args:
            index: The index whose Axis is of interest.

        Returns:
            The ith Axis in the rotation order.

        Raises:
            ValueError: If the index is not in the range (0, 2).
        """
        if index < 0 or index > 2:
            raise ValueError("Index value out of range.")
        return self._rotation_order[index]

    def __str__(self):
        """
        Returns a string representation of the rotation order.

        Returns:
            A string with the axes names of the rotation order.
        """
        return str([str(i) for i in self._rotation_order])

    def __iter__(self):
        self._n = 0
        return self

    def __next__(self):
        if self._n <= 2:
            axis = self._rotation_order[self._n]
            self._n += 1
            return axis
        else:
            raise StopIteration


# Euler rotation orders

XYZ = EulerOrder(Axis.X_AXIS, Axis.Y_AXIS, Axis.Z_AXIS)
XZY = EulerOrder(Axis.X_AXIS, Axis.Z_AXIS, Axis.Y_AXIS)
YXZ = EulerOrder(Axis.Y_AXIS, Axis.X_AXIS, Axis.Z_AXIS)
YZX = EulerOrder(Axis.Y_AXIS, Axis.Z_AXIS, Axis.X_AXIS)
ZXY = EulerOrder(Axis.Z_AXIS, Axis.X_AXIS, Axis.Y_AXIS)
ZYX = EulerOrder(Axis.Z_AXIS, Axis.Y_AXIS, Axis.X_AXIS)
XYX = EulerOrder(Axis.X_AXIS, Axis.Y_AXIS, Axis.X_AXIS)
XZX = EulerOrder(Axis.X_AXIS, Axis.Z_AXIS, Axis.X_AXIS)
YXY = EulerOrder(Axis.Y_AXIS, Axis.X_AXIS, Axis.Y_AXIS)
YZY = EulerOrder(Axis.Y_AXIS, Axis.Z_AXIS, Axis.Y_AXIS)
ZXZ = EulerOrder(Axis.Z_AXIS, Axis.X_AXIS, Axis.Z_AXIS)
ZYZ = EulerOrder(Axis.Z_AXIS, Axis.Y_AXIS, Axis.Z_AXIS)
