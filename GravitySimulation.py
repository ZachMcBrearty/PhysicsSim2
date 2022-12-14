from itertools import combinations
from time import time
import numpy as np
from io import TextIOWrapper as _TextIOWrapper
from random import random
import matplotlib.pyplot as plt

from Animate import graphR, graphRMSr, setAnimate, animateFile, graphEnergies


#G = 1
system_ = "Solar"

if system_ == "Solar":
    G = 6.67 * 10**-11 # m^3.s^-2.kg^-1 = c^2.s.kg^-1

    M = 1.988 * 10**30 # kg, Sun

    m_earth = 5.9723 * 10**24 # kg
    aph_ear = 1.5210 * 10**11 # m
    sem_ear = 1.4960 * 10**11 # m

    solarscale = aph_ear # 1 A.U.

    smooth = 0.001

    scale = (sem_ear, "A.U.")
    timescale = (365.25*24*60*60, "years")
    DEFAULTFILE = "SolarCollapse.bin"
    p=101
    fs=365
    tracelength = 10
else:
    G = 1
    globscale = 20 # r = scale/2 = 1
    scale = (globscale, "arb unit")
    timescale = (1, "s")
    smooth = 0.05

    DEFAULTFILE = "GlobSim.bin"
    p=100
    fs=1
    tracelength = 5
    
timestepB = 0.01
logscale = 1 / np.log(scale[0]/timestepB + 1)

eps = smooth * scale[0]

def calcMasslessForce(pair):
    A, B = pair
    d = A[0] - B[0]
    dabs = np.sqrt(np.sum(d**2))
    F = (G / (dabs**2 + eps**2)**(3/2)) * d # bracket so one number is multiplied onto the array
    return F, dabs

def vCalcMasslessForce(pairs):
    As = pairs[:, 0, 0]
    Bs = pairs[:, 1, 0]
    ds = As - Bs
    dabss = np.sqrt(np.sum(ds**2, axis=1))
    Fs = (G / (dabss**2 + eps**2)**(3/2))[...,None]*ds
    return Fs, np.min(dabss)

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
        dataset = self.Particles
        pij = np.array(np.triu_indices(len(dataset), 1)).T.reshape(-1, 2)
        pairs = dataset[pij]

        Fs, self.mindist = vCalcMasslessForce(pairs)
        # combinations of i and j to give pairs of numbers.
        for n, (i, j) in enumerate(pij):
            self.Particles[i][2] -= Fs[n] * self.ParticleMasses[j] # masses multiplied here to avoid
            self.Particles[j][2] += Fs[n] * self.ParticleMasses[i] # multiplying then dividing
        self.CurTimeStep = self.MaxTimeStep * (np.log(self.minDist/timestepB + 1) * logscale + .001)

    def Update(self) -> None:
        """Update the positions and velocities of the particles in the system based
        on the current time step"""  
        updateMatrix = np.array([[1, self.CurTimeStep, self.CurTimeStep**2],[0, 1, self.CurTimeStep],[0, 0, 0]])
        self.Particles = updateMatrix @ self.Particles
        self.minDist = scale[0]
        self.time += self.CurTimeStep

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
    def Leapfrog(self):
        # offset velocity by -1/2 dt
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
        for p in self.Particles:
            p[1] -= p[2] * self.CurTimeStep / 2

vel = lambda aph, sem: (G * M * (2/aph - 1/sem))**.5
def genRandomPosVel():
    
    r_aph = scale[0] * (0.9+0.2*random())
    r_sem = r_aph #* (random()/4 + 0.75)
    
    theta = (0.45+0.1*random()) * np.pi
    phi = random() * 2 * np.pi
    x = r_aph * np.sin(theta) * np.cos(phi)
    y = r_aph * np.sin(theta) * np.sin(phi)
    z = r_aph * np.cos(theta)
    V_visVisa = vel(r_aph, r_sem)
    theta = np.pi/2
    phi = (phi + np.pi/2) % (2*np.pi)
    vx = V_visVisa * np.cos(phi)
    vy = V_visVisa * np.sin(phi)
    vz = 0
    return x, y, z, vx, vy, vz

def solarCollapse(n=100):
    a = System(0.1 * 60 * 60, file=DEFAULTFILE) # 0.1 hour, in s
    a.AddParticle(M, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) # Sun
    
    for _ in range(n):
        x, y, z, vx, vy, vz = genRandomPosVel()
        a.AddParticle(m_earth/100, x, y, z, vx, vy, vz) # small mass particle

    a.Record()
    a.Leapfrog()
    a.Update()
    n = 3650
    for t in range(n):
        if 100*(t)/n % 10 == 0:
            print(100*(t)/n,"%", end=" ",flush=True)
        a.Record()
        a.doTimestep(24*60*60)
    print("100%", end=" ", flush=True)
    a.Record()
    a.File.close()

def globClust():
    # globular cluster
    a = System(0.2, file=DEFAULTFILE)
    p = 100
    for _ in range(p):
        # sphere radius scale / 2
        R = scale[0] * (1+0.2*random())
        theta = random() * 2 * np.pi
        phi = random() * np.pi
        x = R * np.sin(theta) * np.cos(phi)
        y = R * np.sin(theta) * np.sin(phi)
        z = R * np.cos(theta)
        a.AddParticle(0.01, x, y, z)
    n=1000
    for t in range(n):
        if 100*(t+1)/n % 10 == 0:
            print(100*(t+1)/n,"%", end=" ")
        a.Record()
        a.doTimestep(0.2)
    a.File.close()

def randomP():
    x,y,z,vx,vy,vz = genRandomPosVel()
    return np.array([[x,y,z], [vx,vy,vz], [0.0,0.0,0.0]])

if __name__=="__main__":
    from time import perf_counter
    # rms radius and median radius -> half mass radius for glob
    t0 = perf_counter()
    solarCollapse(100)
    t1 = perf_counter()
    print(t1 - t0, "s")
    setAnimate(widthheight_=scale[0]*1.2, scale_=scale, 
               timescale_=timescale, tracelength_=tracelength)
    # animateFile(DEFAULTFILE, p=51, frameskip=1, repeat=False, ax=(0,1))
    animateFile(DEFAULTFILE, p=51, frameskip=1, repeat=False, ax=(1,2))
    # animateFile(DEFAULTFILE, p=51, frameskip=1, repeat=False, ax=(0,2))
    