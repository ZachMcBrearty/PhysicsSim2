from numpy import degrees, radians, cos, sin, arccos, arctan2

class Vector():
    def __init__(self, x, y, z=0):
        '''inherently 3d vector, z can be left off to make a 2d vector'''
        self.x = x
        self.y = y
        self.z = z

    def rotateTheta(self, theta):
        theta = radians(theta)
        x = self.x*cos(theta) - self.y*sin(theta)
        y = self.x*sin(theta) + self.y*cos(theta)
        return Vector(x, y)

    def rotatePhi(self, phi):
        phi = radians(phi)
        r = (self.x**2+self.y**2+self.z**2)**0.5
        theta = arccos(self.z / r)
        phi = arctan2(self.y, self.x) + phi
        theta = degrees(theta)
        phi = degrees(phi)
        return Vector.fromspherical(r, theta, phi)

    @classmethod
    def fromspherical(cls, r, theta, phi):
        theta = radians(theta)
        phi = radians(phi)
        return r * cls(cos(theta)*sin(phi),
            sin(theta)*sin(phi),
            cos(phi))

    def __str__(self):
        return f"({self.x}, {self.y}, {self.z})"
    __repr__=__str__

    def __add__(self, other):
        if isinstance(other, Vector):
            x = self.x + other.x
            y = self.y + other.y
            z = self.z + other.z
            return Vector(x, y, z)
        else:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, Vector):
            x = self.x - other.x
            y = self.y - other.y
            z = self.z - other.z
            return Vector(x, y, z)
        else:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            x = self.x * other
            y = self.y * other
            z = self.z * other
            return Vector(x, y, z)
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, (int, float)):
            x = self.x * other
            y = self.y * other
            z = self.z * other
            return Vector(x, y, z)
        else:
            return NotImplemented

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            if other == 0:
                return Vector(0, 0, 0)
            else:
                x = self.x / other
                y = self.y / other
                z = self.z / other
                return Vector(x, y, z)
        else:
            return NotImplemented

    def __rtruediv__(self, other):
        if isinstance(other, (int, float)):
            if self.x == 0:
                return Vector(0, 0, 0)
            else:
                x = other / self.x
                y = other / self.y
                z = other / self.z 
                return Vector(x, y, z)
        else:
            return NotImplemented  

    def __abs__(self):
        return (self.x**2 + self.y**2 + self.z**2)**0.5

    def normalise(self):
        return self / abs(self)