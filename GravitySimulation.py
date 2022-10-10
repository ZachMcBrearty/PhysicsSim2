from glob import glob
from itertools import combinations
import numpy as np
from io import TextIOWrapper as _TextIOWrapper
from random import random

from Animate import animateFile, graphEnergies

G = 6.67 * 10**-11 # m^3.s^-2.kg^-1 = c^2.s.kg^-1
M = 1.988 * 10**30 # kg
m = 5.97 * 10**24 # kg
aph = 1.5210 * 10**11 # m
sem = 1.4960 * 10**11 # m

solarscale = aph
globscale = 1
scale = globscale
timestepB = 100

logscale = 1 / np.log(scale/timestepB + 1)

# dampening of the system
eps = 0.1 * scale

# Default file to use if none is specified
DEFAULTFILE = "GravitySim.bin"

class System:
    def __init__(self, MaxTimeStep:float=100, file=None):
        '''Gravity Simulation system, add particles using .AddParticle,
        MaxTimeStep is the largest step in time taken in the simulation, variable based on distance
        file is the name or open fileobject of the binary datafile to be used,
        if none the "DEFAULTFILE" value is used'''
        self.time = 0
        self.MaxTimeStep = MaxTimeStep
        self.CurTimeStep = MaxTimeStep
        self.ParticleMasses = None
        self.Particles = None
        if file is None: 
            self.File = open(DEFAULTFILE, "wb")
        elif isinstance(file, str):
            self.File = open(file, "wb")
        elif isinstance(file, _TextIOWrapper):
            self.File = file
        self.minDist = scale
        self.rec = np.array([])

    def Interaction(self) -> None:
        """Pair up each particle in self.Particles, find the distance between them
        and calculate and add the acceleration they experience due to gravity"""
        l = len(self.Particles)
        # combinations of i and j to give pairs of numbers.
        for i, j in combinations(range(l), 2):
            A = self.Particles[i]
            B = self.Particles[j]
            d = A[0] - B[0]
            dabs = np.sqrt(np.sum(d**2))
            if self.minDist is None or dabs < self.minDist:
                self.minDist = dabs # minimum distance in the simulation to calculate timestep
            F = (G / (dabs**2 + eps**2)**(3/2)) * d # bracket so one number is multiplied onto the array
            A[2] -= F * self.ParticleMasses[j] # masses multiplied here to avoid
            B[2] += F * self.ParticleMasses[i] # multiplying then dividing
        # + 0.001 to avoid having timestep go completely to 0 (ie when 2 particles are overlapping)
        # important in the 1D case
        self.CurTimeStep = self.MaxTimeStep * (np.log(self.minDist/timestepB + 1) * logscale + .001)

    def Update(self) -> None:
        """Update the positions and velocities of the particles in the system based
        on the current time step"""  
        updateMatrix = np.array([[1, self.CurTimeStep, 0],[0, 1, self.CurTimeStep],[0, 0, 0]])
        # adds velocity to position, acceleration to velocity and then sets acceleration to 0
        self.Particles = updateMatrix @ self.Particles
        # reset the minimum distance
        self.minDist = scale
        self.time += self.CurTimeStep

    def doTimestep(self, tmin:float = None) -> None:
        """Call self.Interaction then self.Update, variable timestep length"""
        if tmin is None:
            self.Interaction()
            self.Update()
        else:
            t = 0
            while t < tmin:
                self.Interaction()
                if t + self.CurTimeStep > tmin:
                    # dont overshoot the desired timestep
                    self.CurTimeStep = tmin - t
                self.Update()
                t += self.CurTimeStep

    def AddParticle(self, m:float, x:float, y:float, z:float, vx=0.0, vy=0.0, vz=0.0):
        """Add a particle with mass m at coordinates (x, y, z) 
        with optional velocity (vx, vy, vz)"""
        # these 2 arrays should always be the same length
        if self.ParticleMasses is None:
            self.ParticleMasses = np.array([m])
        else:
            self.ParticleMasses = np.append(self.ParticleMasses, m)

        if self.Particles is None:
            self.Particles = np.array([[[x,y,z], [vx,vy,vz], [0.0,0.0,0.0]]])
        else:
            self.Particles = np.append(self.Particles, [[[x,y,z], [vx,vy,vz], [0.0,0.0,0.0]]], axis=0)
    
    def AddMultipleParticles(self, particles):
        '''must be in the format [p1, p2, ...]
        where pn = [[x, y, z], [vx, vy, vz], [ax, ay, az]]'''
        if self.Particles is None:
            self.Particles = np.array(particles)
        else:
            self.Particles = np.append(self.Particles, particles)

    def kineticEnergy(self):
        """Find and sum all the kinetic energies of the particles in the system"""
        # another matrix multiplication to save time
        return 0.5 * np.sum(self.ParticleMasses.T @ self.Particles[:,1]**2)

    def gravitationalEnergy(self):
        """Find and sum all the graviational energies of the pairs of particles in the system"""
        l = len(self.Particles)
        Et = 0
        # combinations of i and j to give pairs of numbers.
        for i, j in combinations(range(l), 2):
            A = self.Particles[i]
            B = self.Particles[j]
            d = A[0] - B[0]
            dabs = np.sqrt(np.sum(d**2))
            # Eps included here due to change in force changing energy
            Et -= (G * self.ParticleMasses[i] * self.ParticleMasses[j] / (dabs**2 + eps**2)**(1/2))
        return Et

    def Record(self):
        """Record the state of the system to self.File
        state -> time, energy, position"""
        np.array(self.rec).tofile(self.File)
        self.rec = np.array([])

    def append(self):
        kin = self.kineticEnergy()
        grav = self.gravitationalEnergy()
        self.rec = np.append(self.rec, [self.time, kin, grav] + list(self.Particles[:, 0].flatten()))

def solarExample():
    a = System(0.1 * 60 * 60, file="Solar.bin") # 0.1 hour = 6 min = 360s
    vplan = (G * M * (2/aph - 1/sem))**0.5
    # motion of the sun is to effectively conserve momentum so 
    # the system stays fixed at (0,0)
    a.AddParticle(M, 0.0, 0.0, 0.0, 0.0, -m / M * vplan, 0.0) # Sun
    a.AddParticle(m, aph, 0.0, 0.0, 0.0, vplan, 0.0)
    # #print(a.Particles[:,:,0])
    for q in range(365*4 + 1):
        # record first to get t=0
        a.append()
        if q % 10 == 0:
            a.Record()
        # 1/4 of a day to include the .25 in 365.25 days in a year
        a.doTimestep(6*60*60)
    # Record after to get the final timestep
    a.append()
    a.Record()
    # framesskip = 1 day
    animateFile("Solar.bin", framesskip=4, repeat=False,
            scale_=(solarscale, "A.U."), p=2, timescale_=(365.25*24*60*60, "years"))
    graphEnergies("Solar.bin", change=False, p=2)
    graphEnergies("Solar.bin", change=True,  p=2)
    a.File.close()

def globRand():
    return 2*(random()-0.5)

def globularExample():
    # globular cluster
    p=20
    a = System(10, file="Glob.bin")
   
    for _ in range(p):
        a.AddParticle(100, globRand(), globRand(), globRand())

    for t in range(10000):
        a.append()
        if t % 100 == 0:
            print(100*t/10000,"%")
            a.Record()
        a.doTimestep(10)

    animateFile("Glob.bin", framesskip=10, repeat=False,
            scale_=(globscale, "m"), p=p, timescale_=(1, "s"))
    graphEnergies("Glob.bin", p=p)
    graphEnergies("Glob.bin", True, p=p)
    a.File.close()

if __name__=="__main__":
    # solarExample()
    globularExample()