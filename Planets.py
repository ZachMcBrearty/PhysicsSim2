#! python
from itertools import combinations
from math import comb
from numpy import array, arctan2, arccos, real
from numpy import radians, cos, sin, degrees
from random import randint
from datetime import datetime
import matplotlib.pyplot as plt


from Animate import graphEnergies, animateFile

now = datetime.now()
now = now.strftime(f"%M-%H-%d-%m-%y")
filename = f"PlanetsPositionsEnergy{now}.txt"
f = open(filename, "w")
t = 0
# timestep = 1 * 24 * 60 * 60 # seconds, 1 day
timestep = 10_000_000 * 365.25 * 24 * 60 * 60 # seconds, 1 year
# maxtime = 1 * 365.25 * 24 * 60 * 60 # seconds, 5 year
maxtime = 1_00_000_000_000 * 365.25 * 24 * 60 * 60 # seconds, 100,000 years
stepcount = maxtime / timestep + 1
print(stepcount)
frameskip = 25

# Physical Constants
G = 6.67 * 10**-11 # N m^2 kg^-2
soften = 0.1 # used to ensure close interacting bodies do not accelerate infinitely

globClusterSize = 500 * 10**16 # m
solarSize = 2*10**11 # m
scale = globClusterSize

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

class Planet():
    def __init__(self, name:str , mass: float, radius:float,
     colour:tuple, pos:Vector):
        self.name = name
        self.mass = mass
        self.radius = radius
        self.colour = colour
        self.pos = pos
        self.vel = Vector(0,0,0)
        self.acc = Vector(0,0,0)
    @classmethod
    def fromParent(cls, name, rad, ownmass, parent, colour,
     semimajoraxis, aphelion, progression, inclination):
        pos = Vector.fromspherical(aphelion, progression, inclination)
        pl = cls(name, ownmass, rad, colour, parent.pos + pos)
        
        if semimajoraxis == 0:
            vel = Vector(0,0,0)
        else:
            # the velocity of the orbit is at a right angle to the position vector
            vel = pos.normalise()
            vel = vel.rotateTheta(90) # degrees
            # this equation only works when the parent mass is much larger than child
            #vel *= (G * parent.mass / semimajoraxis)**0.5
            # this equation only works for all non-extreme orbits
            vel *= (G * parent.mass * (2 / aphelion - 1 / semimajoraxis)) ** 0.5
        pl.vel = vel + parent.vel
        return pl

    def update(self):
        self.vel += self.acc * timestep
        self.pos += self.vel * timestep
        self.acc = Vector(0,0,0)

def EarthSunSystem():
    Sun = {"Name": "Sun",
        "Mass": 1.989 * 10**30, # kg
        "Radius": 6.957 * 10**8, # m
        "Colour": (255, 255, 0),
        } # Yellow
    Earth = {"Name": "Earth",
            "Mass": 5.9723 * 10**24, # kg
            "Radius": 6.371 * 10**6, # m
            "Colour": (0, 0, 255),
            "Semimajor": 149.60 * 10**9,# m, average distance in orbit
            "Aphelion": 152.10 * 10**9,# m, highest point in orbit
            "Inclination": 90 # degrees
            } # Blue
    system = [Planet(Sun["Name"], Sun["Mass"],Sun["Radius"], Sun["Colour"], Vector(0,0,0))]
    system.append(Planet.fromParent(Earth["Name"],
    Earth["Radius"], Earth["Mass"], system[0], 
    Earth["Colour"], Earth["Semimajor"], Earth["Aphelion"],
    randint(0,360), Earth["Inclination"]))
    return system

def randomSystem2D(N, range_=10**9):
    return [Planet(str(q), 10**24, 10**6, (0,0,0), 
            Vector(randint(-range_, range_), randint(-range_, range_), 0)) for q in range(N)]

def randomSystem3D(N, range_=10**9):
    return [Planet(str(q), 10**24, 10**6, (0,0,0), 
            Vector(randint(-range_, range_), randint(-range_, range_), randint(-range_, range_))) for q in range(N)]

def findForce(A, B):
    d = A.pos - B.pos
    return d * G * A.mass * B.mass / (abs(d)**2 + (10**6)**2)**(3/2)

def interact(system):
    pairs = combinations(system, 2)
    for A, B in pairs:
        F = findForce(A, B)
        A.acc -= F / A.mass
        B.acc += F / B.mass

def update(system):
    for s in system:
        s.update()

def calcEnergy(system):
    kinetic = 0
    gravity = 0
    
    for q in system:
        kinetic += 1/2 * q.mass * abs(q.vel)**2

    for a, b in combinations(system, 2):
        d = a.pos - b.pos
        gravity -= G * a.mass * b.mass / (abs(d)**2 + (10**10)**2)**0.5

    return str(kinetic), str(gravity), str(kinetic+gravity)

def record(system, file):
    # KinEnergy, GravEn, TotalEn # (x1, y1, z1), (x2, y2, z2)...
    energy = " , ".join(calcEnergy(system))
    positions = " * ".join(map(lambda x: str(x.pos), system))
    file.write(str(t)+ " # " + energy + " # " + positions + "\n")

if __name__=="__main__":
    #system = EarthSunSystem()
    system = randomSystem2D(50, scale)
    for step in range(int(stepcount)):
        interact(system)
        update(system)
        record(system, f)
        t += timestep
    f.close()

    animateFile(filename, framesskip=frameskip, repeat=False, scale_=scale*10)
    graphEnergies(filename, True)
    
