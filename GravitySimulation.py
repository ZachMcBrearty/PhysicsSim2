from itertools import combinations
from time import time
import numpy as np
from io import TextIOWrapper as _TextIOWrapper
from random import random
import matplotlib.pyplot as plt

from Animate import setAnimate, animateFile, graphEnergies


#G = 1
system_ = "GlobClust"

if system_ in ["Solar", "SolarPlus", "GlobClust"]:
    if system_ in ["Solar", "SolarPlus"]:
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
            p=2
            fs=365
            tracelength = -1
        else:
            scale = (sem_ear, "A.U.")
            timescale = (365.25*24*60*60, "years")
            DEFAULTFILE = "SolarPlus.bin"
            p=3
            fs=5
            tracelength = -1
    else:
        G = 1
        globscale = 20 # r = scale/2 = 1
        scale = (globscale, "arb unit")
        timescale = (1, "s")
        smooth = 0.05

        DEFAULTFILE = "GlobSim.bin"
        p=20
        fs=1
        tracelength = 5
    
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
            d1 = (A[0] - B[0]) % scale[0]
            d2 = (B[0] - A[0]) % scale[0]
            d = np.array([max(a, b) for a, b in zip(d1, d2)])
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
        self.Particles[:, 0] %= scale[0]
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
            d1 = (A[0] - B[0]) % scale[0]
            d2 = (B[0] - A[0]) % scale[0]
            d = np.array([max(a, b) for a, b in zip(d1, d2)])
            dabs = np.sqrt(np.sum(d**2))
            Et -= (G * self.ParticleMasses[i] * self.ParticleMasses[j] / (dabs**2 + eps**2)**(1/2))
        return Et

    def Record(self):
        """Record the state of the system to self.File
        state -> time, energy, position"""
        kin = self.kineticEnergy()
        grav = self.gravitationalEnergy()
        a = [self.time, kin, grav] + list((self.Particles[:, 0]-scale[0]/2).flatten())
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
        

def solarExample():
    a = System(0.1 * 60 * 60, file=DEFAULTFILE) # 0.1 hour, in s
    vplan = (G * M * (2/aph_ear - 1/sem_ear))**0.5
    # motion of the sun is to effectively conserve momentum so 
    # the system stays fixed at (0,0)
    a.AddParticle(M, 0.0, 0.0, 0.0, 0.0, -m_earth / M * vplan, 0.0) # Sun
    a.AddParticle(m_earth, aph_ear, 0.0, 0.0, 0.0, vplan, 0.0)
    # #print(a.Particles[:,:,0])
    a.Record()
    a.Leapfrog()
    a.Update()
    n = 36525
    for t in range(n):
        if 100*(t+1)/n % 10 == 0:
            print(100*(t+1)/n,"%", end=" ")
        a.Record()
        a.doTimestep(24*60*60)
    a.Record()
    a.File.close()

def variableTimestep():
    percentchange = []
    percentchangeLeap = []
    # 1 day, 6 hours, 1 hour, 10 minutes, 1 minute
    ts = [24*60*60, 6*60*60, 60*60, 10*60, 60]
    tmax = 365*24*60*60
    vplan = (G * M * (2/aph_ear - 1/sem_ear))**0.5
    for t in ts:
        print(f"T: {t}", end=" ")
        a = System(t, file=DEFAULTFILE)
        a.AddParticle(M, 0.0, 0.0, 0.0, 0.0, -m_earth / M * vplan, 0.0) # Sun
        a.AddParticle(m_earth, aph_ear, 0.0, 0.0, 0.0, vplan, 0.0)  
        e0 = a.gravitationalEnergy() + a.kineticEnergy()
        a.doTimestep(tmax)
        ef = a.gravitationalEnergy() + a.kineticEnergy()
        percentchange.append(100*abs((ef - e0) / e0))

        b = System(t, file=DEFAULTFILE)
        b.AddParticle(M, 0.0, 0.0, 0.0, 0.0, -m_earth / M * vplan, 0.0) # Sun
        b.AddParticle(m_earth, aph_ear, 0.0, 0.0, 0.0, vplan, 0.0)  
        e0 = b.gravitationalEnergy() + b.kineticEnergy()
        b.Leapfrog()
        b.Update()
        b.doTimestep(tmax)
        ef = b.gravitationalEnergy() + b.kineticEnergy()
        percentchangeLeap.append(100*abs((ef - e0) / e0))

        a.File.close()
        b.File.close()
    print()
    plt.loglog(ts, percentchange, label="No Leapfrog")
    plt.loglog(ts, percentchangeLeap, label="Leapfrog")
    plt.legend()
    plt.show()

def globClust():
    # globular cluster
    a = System(0.1, file=DEFAULTFILE)
    for _ in range(p):
        # sphere radius scale / 2
        x = scale[0] * random()
        y = scale[0] * random()
        z = scale[0] * random()
        a.AddParticle(0.01, x, y, z)

    n=1000
    for t in range(n):
        if 100*(t+1)/n % 10 == 0:
            print(100*(t+1)/n,"%", end=" ")
        a.Record()
        a.doTimestep(1)
    a.File.close()

def earthSunJupiter():
    #DEFAULTFILE = "SolarPlus20yr.bin"
    a = System(1 * 60 * 60, file=DEFAULTFILE) # 1 hour, in s
    v_earth = (G * M * (2/aph_ear - 1/sem_ear))**0.5
    v_jup = (G * M * (2/aph_jup - 1/sem_jup))**0.5
    # motion of the sun is to effectively conserve momentum so 
    # the system stays fixed at (0,0)
    # SUN
    a.AddParticle(M, 0.0, 0.0, 0.0, 0.0, -m_earth / M * v_earth + m_jup / M * v_jup, 0.0) # Sun
    # EARTH
    a.AddParticle(m_earth, aph_ear, 0.0, 0.0, 0.0, v_earth, 0.0)
    # JUPTIER
    a.AddParticle(m_jup, -aph_jup, 0, 0, 0, -v_jup, 0)
    total = int(365.25 * 400) #int(4333 * 100)
    a.Record()
    a.Leapfrog()
    a.Update()
    for q in range(total):
        if (100*((q+1) / total)) % 10 == 0:
            print(100*(q+1) / total, "%", end=" ")
        a.Record()
        a.doTimestep(24*60*60)
    a.Record()
    # framesskip = 20 days
    a.File.close()  

if __name__=="__main__":
    if True:
        t0 = time()
        if system_ == "Solar":
            solarExample()
        elif system_ == "SolarPlus":
            earthSunJupiter()
        else:
            globClust()
        t1 = time()
        te = (t1-t0)
        if te > 60*60:
            print("Time taken:", te/3600, "hours")
        elif te > 60:
            print("Time taken:", te/60, "seconds")
        else:
            print("Time taken:", te, "seconds")
        

    setAnimate(widthheight_=scale[0]*1.2, scale_=scale, 
               timescale_=timescale, tracelength_=tracelength)
    animateFile(DEFAULTFILE, p=p, 
                frameskip=fs, repeat=True)
    graphEnergies(DEFAULTFILE, False, p=p)
    graphEnergies(DEFAULTFILE, True,  p=p)

    # variableTimestep()