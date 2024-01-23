from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any


@dataclass(frozen=True)
class Point:
    x: Any
    y: Any


class AnalyticFunction(ABC):

    @abstractmethod
    def compute(self, *args, **kwargs):
        pass

    def computePoint(self, *args, **kwargs):
        """The first positional argument will be taken to be the x-value for the Point."""

        xValue = args[0]
        yValue = self.compute(*args, **kwargs)

        return Point(xValue, yValue)


@dataclass
class Boundary:
    lower: field(init=False)
    upper: field(init=False)
    function: AnalyticFunction

    def __init__(self, point1: Point, point2: Point, function: AnalyticFunction):
        if point1.x < point2.x:
            self.lower = point1
            self.upper = point2
        else:
            self.upper = point1
            self.lower = point2
        self.function = function

    def bifurcate(self, *args, **kwargs):
        """Bifurcate and return two new Boundary objects with the middle Point as the
        upper and lower point for the left and right boundaries respectively. Arguments
        are passed to the AnalyticFunction.computePoint method and should be passed exactly as
        they would be to that method."""

        middleValue = self.function.computePoint(*args, **kwargs)
        leftBoundary = Boundary(self.lower, middleValue, self.function)
        rightBoundary = Boundary(middleValue, self.upper, self.function)

        return leftBoundary, rightBoundary

    def range(self):
        """Returns the difference between the upper and lower x values. This value
        is always positive."""

        return self.upper.x - self.lower.x

    def difference(self):
        """Returns the absolute value of the y values. This value is always positive."""

        return abs(self.upper.y - self.lower.y)

    def hasSameSign(self):
        """Returns True if the signs of the y values are the same. If they are different
        or at least one of the values is 0, returns False."""

        return self.lower.y * self.upper.y > 0
