from enum import Enum


# enum for general axis directions, empty_axis has been useful before (2 rotations)
class Axis(Enum):
    X_AXIS = 0
    Y_AXIS = 1
    Z_AXIS = 2
    EMPTY_AXIS = 3


class EulerOrder:

    def __init__(self, axis1: Axis, axis2: Axis, axis3: Axis):
        self._rotation_order = (axis1, axis2, axis3)
        self._first_rotation = axis1
        self._second_rotation = axis2
        self._third_rotation = axis3

    def __getitem__(self, i: int) -> Axis:
        if i < 0 or i > 2:
            raise ValueError("Index value out of range.")
        return self._rotation_order[i]

    def __str__(self):
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


class Order:
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
    X = EulerOrder(Axis.X_AXIS, Axis.EMPTY_AXIS, Axis.EMPTY_AXIS)
    Y = EulerOrder(Axis.Y_AXIS, Axis.EMPTY_AXIS, Axis.EMPTY_AXIS)
    Z = EulerOrder(Axis.Z_AXIS, Axis.EMPTY_AXIS, Axis.EMPTY_AXIS)
