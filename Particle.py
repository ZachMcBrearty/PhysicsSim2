from Vector import Vector 
from Constants import *
from random import randint

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
        self.local = timestep
    def setLocal(self, lts):
        if lts < self.local:
            self.local = lts
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
        self.pos += self.vel * self.local
        self.vel += self.acc * self.local
        self.acc = Vector(0,0,0)

        

