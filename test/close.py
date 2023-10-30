from math import isclose

from pyevspace import Vector


def vectorIsClose(lhs: Vector, rhs: Vector, places: int = 9) -> bool:
    epsilon = pow(10, -places)
    for left, right in zip(lhs, rhs):
        if abs(left - right) > epsilon:
            return False

        # if isclose(left, right, rel_tol=rel_tol, abs_tol=abs_tol) is False:
        #     return False

    return True
