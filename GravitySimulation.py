from itertools import combinations
from time import time
import numpy as np
from io import TextIOWrapper as _TextIOWrapper
from random import random

from Animate import setAnimate, animateFile, graphEnergies

#G = 1
system_ = "GlobClust"

if system_ in ["Solar", "SolarPlus", "GlobClust"]:
    if system_ in ["Solar", "SolarPlus"]:
        c = 3 * 10**8 # m.s^-1
        G = 6.67 * 10**-11 # m^3.s^-2.kg^-1 = c^2.s.kg^-1
    
        M = 1.988 * 10**30 # kg, Sun

        m_earth = 5.9723 * 10**24 # kg
        aph_ear = 1.5210 * 10**11 # m
        sem_ear = 1.4960 * 10**11 # m

        m_jup = 1.8986 * 10**27 # kg
        aph_jup = 816.62 * 10**9 # m
        sem_jup = 778.57 * 10**9 # m

        solarscale = aph_ear # 1 A.U.
        solarplusscale = aph_jup

        smooth = 0

        if system_ == "Solar":
            scale = (sem_ear, "A.U.")
            timescale = (365.25*24*60*60, "years")
            DEFAULTFILE = "Solar.bin"
        else:
            scale = (sem_ear, "A.U.")
            timescale = (365.25*24*60*60, "years")
            DEFAULTFILE = "SolarPlus.bin"
    else:
        G = 1
        globscale = 2
        scale = (globscale, "arb unit")
        timescale = (1, "s")
        smooth = 0.1

        DEFAULTFILE = "GlobSim.bin"
    
timestepB = 100
logscale = 1 / np.log(scale[0]/timestepB + 1)

eps = smooth * scale[0]

class System:
    def __init__(self, MaxTimeStep=100, file=None):
        '''Gravity Simulation system, add particles using .AddParticle,
        MaxTimeStep is the largest step in time taken in the simulation, variable based on distance
        file is the name or fileobject of the binary datafile to be used,
        if none used the "DEFAULTFILE" value'''
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
        self.minDist = scale[0]

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
        self.CurTimeStep = self.MaxTimeStep * (np.log(self.minDist/timestepB + 1) * logscale + .001)

    def Update(self) -> None:
        """Update the positions and velocities of the particles in the system based
        on the current time step"""  
        updateMatrix = np.array([[1, self.CurTimeStep, 0],[0, 1, self.CurTimeStep],[0, 0, 0]])
        self.Particles = updateMatrix @ self.Particles
        self.minDist = scale[0]
        self.time += self.CurTimeStep
        # for particle in self.Particles:
        #     particle = updateMatrix @ particle

    def doTimestep(self, tmin:float = None) -> None:
        """Call self.Interaction then self.Update, variable timestep length"""#
        if tmin is None:
            self.Interaction()
            self.Update()
        else:
            t = 0
            while t < tmin:
                self.Interaction()
                if t + self.CurTimeStep > tmin:
                    self.CurTimeStep = tmin - t
                self.Update()
                t += self.CurTimeStep

    def AddParticle(self, m:float, x:float, y:float, z:float, vx=0.0, vy=0.0, vz=0.0):
        """Add a particle with mass m at coordinates (x, y, z) with optional velocity (vx, vy, vz)"""
        # idk which one works
        # self.Particles.append(np.array([[x,y,z], [vx,vy,vz],[0,0,0]]))
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
            Et -= (G * self.ParticleMasses[i] * self.ParticleMasses[j] / (dabs**2 + eps**2)**(1/2))
        return Et

    def Record(self):
        """Record the state of the system to self.File
        state -> time, energy, position"""
        kin = self.kineticEnergy()
        grav = self.gravitationalEnergy()
        a = [self.time, kin, grav] + list(self.Particles[:, 0].flatten())
        np.array(a).tofile(self.File)

    def print(self):
        print(self.ParticleMasses)
        print(self.Particles)

def solarExample():
    a = System(0.1 * 60 * 60, file=DEFAULTFILE) # 0.1 hour, in s
    vplan = (G * M * (2/aph_ear - 1/sem_ear))**0.5
    # motion of the sun is to effectively conserve momentum so 
    # the system stays fixed at (0,0)
    a.AddParticle(M, 0.0, 0.0, 0.0, 0.0, -m_earth / M * vplan, 0.0) # Sun
    a.AddParticle(m_earth, aph_ear, 0.0, 0.0, 0.0, vplan, 0.0)
    # #print(a.Particles[:,:,0])
    for _ in range(366):
        a.Record()
        a.doTimestep(24*60*60)
    # framesskip = 20 days
    setAnimate(scale_=scale, timescale_=timescale)
    animateFile(DEFAULTFILE, p=2, frameskip=1, repeat=False)
    graphEnergies(DEFAULTFILE, False)
    graphEnergies(DEFAULTFILE, True)
    a.File.close()

def globClust():
    # globular cluster
    a = System(0.01, file=DEFAULTFILE)
    p = 100
    for _ in range(p):
        # sphere radius scale / 2
        R = scale[0] * (random()-0.5)
        theta = random() * 2 * np.pi
        phi = random() * np.pi
        x = R * np.sin(theta) * np.cos(phi)
        y = R * np.sin(theta) * np.sin(phi)
        z = R * np.cos(theta)
        a.AddParticle(0.01, x, y, z)
    n=100
    for t in range(n):
        if 100*(t+1)/n % 10 == 0:
            print(100*(t+1)/n,"%")
        a.Record()
        a.doTimestep(0.1)
    setAnimate(widthheight_=scale[0]*1.2, scale_=scale, 
               timescale_=timescale)
    animateFile(DEFAULTFILE, p=p, 
                frameskip=1, repeat=False)
    graphEnergies(DEFAULTFILE, False, p=p)
    graphEnergies(DEFAULTFILE, True,  p=p)

def earthSunJupiter():
    DEFAULTFILE = "SolarPlus20yr.bin"
    # a = System(1 * 60 * 60, file=DEFAULTFILE) # 1 hour, in s
    # v_earth = (G * M * (2/aph_ear - 1/sem_ear))**0.5
    # v_jup = (G * M * (2/aph_jup - 1/sem_jup))**0.5
    # # motion of the sun is to effectively conserve momentum so 
    # # the system stays fixed at (0,0)
    # # SUN
    # a.AddParticle(M, 0.0, 0.0, 0.0, 0.0, -m_earth / M * v_earth + m_jup / M * v_jup, 0.0) # Sun
    # # EARTH
    # a.AddParticle(m_earth, aph_ear, 0.0, 0.0, 0.0, v_earth, 0.0)
    # # JUPTIER
    # a.AddParticle(m_jup, -aph_jup, 0, 0, 0, -v_jup, 0)
    # total = int(365.25 * 100) #int(4333 * 100)
    # for q in range(total):
    #     if (100*((q+1) / total)) % 10 == 0:
    #         print(100*(q+1) / total, "%")
    #     a.Record()
    #     a.doTimestep(24*60*60)
    # a.Record()
    # framesskip = 20 days
    setAnimate(solarplusscale*1.2, scale, timescale)
    animateFile(DEFAULTFILE, p=3, frameskip=365, repeat=False)
    graphEnergies(DEFAULTFILE, False, p=3)
    graphEnergies(DEFAULTFILE, True,  p=3)
    # a.File.close()  

if __name__=="__main__":
    # earthSunJupiter()
    globClust()