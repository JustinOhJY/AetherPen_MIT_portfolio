import math

class Quaternion:
    def __init__(self, s: float, i: float, j: float, k: float):
        self._s = s
        self._i = i
        self._j = j
        self._k = k

    @property
    def s(self) -> float:
        return self._s

    @s.setter
    def s(self, s: float):
        self._s = s

    @property
    def i(self) -> float:
        return self._i

    @i.setter
    def i(self, i: float):
        self._i = i

    @property
    def j(self) -> float:
        return self._j

    @j.setter
    def j(self, j: float):
        self._j = j

    @property
    def k(self) -> float:
        return self._k

    @k.setter
    def k(self, k: float):
        self._k = k

    @property
    def scaler(self) -> float:
        # Computes the magnitude of the vector (i, j, k)
        return math.sqrt(self._i ** 2 + self._j ** 2 + self._k ** 2)

    @property
    def angle(self) -> float:
        # Clamp the value to the valid range for acos
        value = max(-1, min(1, self._s / self.scaler)) if self.scaler != 0 else 0
        return math.degrees(math.acos(value))

    @property
    def s_unit(self) -> float:
        return self._s / self.scaler if self.scaler != 0 else 0

    @property
    def i_unit(self) -> float:
        return self._i / self.scaler if self.scaler != 0 else 0

    @property
    def j_unit(self) -> float:
        return self._j / self.scaler if self.scaler != 0 else 0

    @property
    def k_unit(self) -> float:
        return self._k / self.scaler if self.scaler != 0 else 0

    @property
    def cijks(self):
        # Returns a list containing the angle and normalized components
        if math.sin(math.radians(self.angle)) == 0:
            return [self.angle, 0, 0, 0]
        return [
            self.angle,
            self.i_unit / math.sin(math.radians(self.angle)),
            self.j_unit / math.sin(math.radians(self.angle)),
            self.k_unit / math.sin(math.radians(self.angle))
        ]

    def __mul__(self, q: "Quaternion") -> "Quaternion":
        result = Quaternion(0, 1, 0, 0)
        result.s = self._s * q._s - self._i * q._i - self._j * q._j - self._k * q._k
        result.i = self._s * q._i + self._i * q._s + self._j * q._k - self._k * q._j
        result.j = self._s * q._j - self._i * q._k + self._j * q._s + self._k * q._i
        result.k = self._s * q._k + self._i * q._j - self._j * q._i + self._k * q._s
        return result

    def inverse(self) -> "Quaternion":
        inverse = Quaternion(self._s, -self._i, -self._j, -self._k)
        return inverse

    def __str__(self) -> str:
        return (f'({self.s}, {self.i}, {self.j}, {self.k}) '
                f'({self.scaler} * (cos {self.angle} + sin {self.angle} * ({self.i}i + {self.j}j + {self.k}k)))')

    def rotate(self, i: float, j: float, k: float, angle: float) -> None:
        # Rotate the quaternion by creating a new rotation quaternion from the given axis and angle
        rad = math.radians(angle)
        q1 = Quaternion(math.cos(rad), i * math.sin(rad), j * math.sin(rad), k * math.sin(rad))
        rotated = q1 * self * q1.inverse()
        self._s = rotated.s
        self._i = rotated.i
        self._j = rotated.j
        self._k = rotated.k


def quaternion_rotation(w: float, x: float, y: float, z: float, v: list) -> None:
    # Rotates vector v in place using the quaternion components.
    new_v = [0, 0, 0]
    new_v[0] = (1 - 2 * y * y - 2 * z * z) * v[0] + (2 * x * y - 2 * w * z) * v[1] + (2 * x * z + 2 * w * y) * v[2]
    new_v[1] = (2 * x * y + 2 * w * z) * v[0] + (1 - 2 * x * x - 2 * z * z) * v[1] + (2 * y * z - 2 * w * x) * v[2]
    new_v[2] = (2 * x * z - 2 * w * y) * v[0] + (2 * y * z + 2 * w * x) * v[1] + (1 - 2 * x * x - 2 * y * y) * v[2]
    v[:] = new_v

def quaternion_rotation1(w: float, x: float, y: float, z: float, v: list) -> None:
    # Rotates vector v in place using the quaternion components.
    new_v = [0, 0, 0]
    new_v[0] = (1 - 2 * y * y - 2 * z * z) * v[0] + (2 * x * y - 2 * w * z) * v[1] + (2 * x * z + 2 * w * y) * v[2]
    new_v[1] = (2 * x * y + 2 * w * z) * v[0] + (1 - 2 * x * x - 2 * z * z) * v[1] + (2 * y * z - 2 * w * x) * v[2]
    new_v[2] = (2 * x * z - 2 * w * y) * v[0] + (2 * y * z + 2 * w * x) * v[1] + (1 - 2 * x * x - 2 * y * y) * v[2]
    v[:] = new_v
