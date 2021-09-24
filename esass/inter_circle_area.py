"""Intersection area of two circles"""

import math
from dataclasses import dataclass
from typing import Tuple


@dataclass
class Circle:
    x: float
    y: float
    r: float

    @property
    def coord(self):
        return self.x, self.y


def find_intersection(c1: Circle, c2: Circle) -> float:
    """Finds intersection area of two circles.

    Returns intersection area of two circles otherwise 0
    """

    d = math.dist(c1.coord, c2.coord)
    rad1sqr = c1.r ** 2
    rad2sqr = c2.r ** 2

    if d == 0:
        # the circle centers are the same
        return math.pi * min(c1.r, c2.r) ** 2

    angle1 = (rad1sqr + d ** 2 - rad2sqr) / (2 * c1.r * d)
    angle2 = (rad2sqr + d ** 2 - rad1sqr) / (2 * c2.r * d)

    # check if the circles are overlapping
    if (-1 <= angle1 < 1) or (-1 <= angle2 < 1):
        theta1 = math.acos(angle1) * 2
        theta2 = math.acos(angle2) * 2

        area1 = (0.5 * theta2 * rad2sqr) - (0.5 * rad2sqr * math.sin(theta2))
        area2 = (0.5 * theta1 * rad1sqr) - (0.5 * rad1sqr * math.sin(theta1))

        return area1 + area2
    elif angle1 < -1 or angle2 < -1:
        # Smaller circle is completely inside the largest circle.
        # Intersection area will be area of smaller circle
        # return area(c1_r), area(c2_r)
        return math.pi * min(c1.r, c2.r) ** 2
    return 0


if __name__ == "__main__":

    @dataclass
    class Test:
        data: Tuple[Circle, Circle]
        expected: float

    tests = [
        Test((Circle(2, 4, 2), Circle(3, 9, 3)), 0),
        Test((Circle(0, 0, 2), Circle(-1, 1, 2)), 7.0297),
        Test((Circle(1, 3, 2), Circle(1, 3, 2.19)), 12.5664),
        Test((Circle(0, 0, 2), Circle(-1, 0, 2)), 8.6084),
        Test((Circle(4, 3, 2), Circle(2.5, 3.5, 1.4)), 3.7536),
        Test((Circle(3, 3, 3), Circle(2, 2, 1)), 3.1416)
    ]

    for test in tests:
        result = find_intersection(*test.data)
        assert math.isclose(result, test.expected, rel_tol=1e-4), f"{test=}, {result=}"
    print("PASSED!!!")