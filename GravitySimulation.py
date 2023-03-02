from itertools import combinations
from time import time
import numpy as np
import copy
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

    m_jup = 1.8986 * 10**27 # kg
    aph_jup = 816.62 * 10**9 # m
    sem_jup = 778.57 * 10**9 # m

    solarscale = aph_ear # 1 A.U.

    smooth = 0.001 # 10^-3

    scale = (sem_ear, "A.U.")
    timescale = (365.25*24*60*60, "years")
    DEFAULTFILE = "SolarCollapseJ1.bin"
    p=102
    fs=365
    tracelength = 2

    frameskip = 1
else:
    G = 1
    globscale = 20 # r = scale/2 = 11
    scale = (globscale, "arb unit")
    timescale = (1, "s")
    smooth = 0.05

    DEFAULTFILE = "GlobSim.bin"
    p=100
    fs=1
    tracelength = 5
    frameskip = 1

eps = smooth * scale[0]
rho = 5000 # kg m^-3
gamma = (3 / (4 * np.pi * rho)) / 4 / smooth # tweak to give same as eps

def calcMasslessForce(pair):
    A, B = pair
    d = A[0] - B[0]
    dabs = np.sqrt(np.sum(d**2))
    F = (G / (dabs**2 + eps**2)**(3/2)) * d # bracket so one number is multiplied onto the array
    return F, dabs

def vCalcMasslessForce(pairs):
    #pairs = [[[[x,y,z],[vx,vy,vz],[ax,ay,az]], 
    #          [[x,y,z],[vx,vy,vz],[ax,ay,az]]],
    #         [..., ...], 
    #         ...]
    # Every pair, first of the pair, position values
    As = pairs[:, 0, 0]
    # Every pair, second of the pair, position values
    Bs = pairs[:, 1, 0]
    ds = As - Bs
    dabss = np.sqrt(np.sum(ds**2, axis=1))
    # let eps prop to radius prop to mass**(1/3)
    # requires mass-full force to be calculated?
    Fs = (G / (dabss**2 + eps**2)**(3/2))[...,None]*ds
    return Fs, dabss

class System:
    def __init__(self, MaxTimeStep=100, file=None):
        '''Gravity Simulation system, add particles using .AddParticle,
        MaxTimeStep is the largest step in time taken in the simulation, variable based on distance
        file is the name or fileobject of the binary datafile to be used,
        if none used the "DEFAULTFILE" value'''
        self.time = 0
        self.MaxTimeStep = MaxTimeStep
        self.updateMatrix = np.array([[1, self.MaxTimeStep, self.MaxTimeStep**2],[0, 1, self.MaxTimeStep],[0, 0, 0]])
        self.ParticleMasses = None
        self.Particles = None
        self.coupled = None
        if file is None: 
            self.File = open(DEFAULTFILE, "wb")
        elif isinstance(file, str):
            self.File = open(file, "wb")
        elif isinstance(file, _TextIOWrapper):
            self.File = file
        self.minDist = scale[0]

    def upperTriIter(self, l):
        q = np.array(np.triu_indices(l, 1)).T.reshape(-1, 2)
        l = None
        for p in q:
            if self.coupled[p[0]] and self.coupled[p[1]]:
                if l is None:
                    l = np.array([p])
                else:
                    l = np.append(l, [p], axis=0)
        return l

    def Interaction(self) -> None:
        """Pair up each particle in self.Particles, find the distance between them
        and calculate and add the acceleration they experience due to gravity"""
        dataset = self.Particles[self.coupled]
        if len(dataset) == 1:
            return
        pij = self.upperTriIter(len(self.Particles))
        pairs = self.Particles[pij]

        Fs, dists = vCalcMasslessForce(pairs)
        coup = copy.deepcopy(self.coupled)
       
        epss = gamma * np.sum(self.ParticleMasses[pij] ** (1/3), axis=1)

        for n, (i, j) in enumerate(pij):
            if not coup[i] or not coup[j]:
                print("NOT WORKING 1")
                continue
            if dists[n] < 2*epss[n]:
                if self.ParticleMasses[i] >= self.ParticleMasses[j]:
                    self.Couple(i, j)
                else:
                    self.Couple(j, i)

        coup = copy.deepcopy(self.coupled)

        ### combinations of i and j to give pairs of numbers.
        for n, (i, j) in enumerate(pij):
            if not self.coupled[i] or not self.coupled[j]:
                # coupled in the previous loop, therefore skip
                # print("NOT WORKING 2")
                continue
            self.Particles[i][2] -= Fs[n] * self.ParticleMasses[j] # masses multiplied here to avoid
            self.Particles[j][2] += Fs[n] * self.ParticleMasses[i] # multiplying then dividing

    def Update(self) -> None:
        """Update the positions and velocities of the particles in the system based
        on the current time step"""  
        self.Particles = self.updateMatrix @ self.Particles
        self.time += self.MaxTimeStep

    def Couple(self, i, j):
        """set the coupled condition on j and set the mass to 0
        add the mass of j to i and add the velocities (conserving energy)"""
        # v_new = p1+p2 / (m1+m2)
        # conserves momentum but not energy
        print(f"t={self.time}; Coupled {i} {self.ParticleMasses[i] :.3} and {j} {self.ParticleMasses[j] :.3}, {np.count_nonzero(self.coupled)-1} particles left", flush=True)
        # print(f"Particle {i}: {self.Particles[i]}")
        # print(f"Particle {j}: {self.Particles[j]}")
    
        # set acceleration to 0
        self.Particles[i][2] = [0, 0, 0]
        # move to center of momentum
        self.Particles[i][1] = (self.ParticleMasses[i]*self.Particles[i][1] + self.ParticleMasses[j]*self.Particles[j][1]) / (self.ParticleMasses[i] + self.ParticleMasses[j])
        # move to Centre of Mass
        self.Particles[i][0] = (self.ParticleMasses[i] * self.Particles[i][0] + self.ParticleMasses[j] * self.Particles[j][0]) / (self.ParticleMasses[i] + self.ParticleMasses[j])
        # combine masses
        self.ParticleMasses[i] += self.ParticleMasses[j]
        # move the other particle out and set the mass to 0, and the coupled condition
        self.Particles[j] = np.array([[10*scale[0], 10*scale[0], 10*scale[0]], [0, 0, 0], [0, 0, 0]])
        self.ParticleMasses[j] = 0
        self.coupled[j] = False
        
    def doTimestep(self, tmin:float = None) -> None:
        """Call self.Interaction then self.Update, variable timestep length"""#
        self.Interaction()
        self.Update()
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
        if self.coupled is None:
            self.coupled = np.array([True])
        else:
            self.coupled = np.append(self.coupled, True)
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
        state -> time, energy, masses, position velocity"""
        kin = self.kineticEnergy()
        grav = self.gravitationalEnergy()
        b = [self.time, kin, grav] + list(self.ParticleMasses) + list(self.Particles[:, 0:2].flatten())
        # a = [self.time, kin, grav] + list(self.Particles[:, 0].flatten())
        # np.array(a).tofile(self.File)
        np.array(b).tofile(self.File)
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
        for p in self.Particles:
            p[1] -= p[2] * self.MaxTimeStep / 2

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

def solarCollapseJup(n=100, File=DEFAULTFILE, numYears=50):
    a = System(10 * 24 * 60 * 60, file=File) # 2 days
    v_jup = vel(aph_jup, sem_jup)
    a.AddParticle(M, 0.0, 0.0, 0.0, 0.0, -m_jup * v_jup / M, 0.0) # Sun
    a.AddParticle(m_jup, aph_jup, 0.0, 0.0, 0.0, +v_jup, 0.0) # Jupiter
    for _ in range(n):
        x, y, z, vx, vy, vz = genRandomPosVel()
        a.AddParticle(m_earth/100, x, y, z, vx, vy, vz) # small mass particle
    a.Leapfrog()
    a.Record()
    a.Update()
    n = int(round(numYears * 36.525, -1))
    for t in range(n):
        if 100*(t)/n % 20 == 0:
            print(100*(t)/n,"%", end=" ",flush=True)
        a.Record()
        a.doTimestep()
    print("100%", end=" ", flush=True)
    a.Record()
    a.File.close()

def solarCollapseNoJup(n=100, File=DEFAULTFILE, numYears=50):
    a = System(10 * 24 * 60 * 60, file=File) # 10 days
    a.AddParticle(M, 0.0, 0.0, 0.0, 0.0, 0, 0.0) # Sun
    for _ in range(n):
        x, y, z, vx, vy, vz = genRandomPosVel()
        a.AddParticle(m_earth/100, x, y, z, vx, vy, vz) # small mass particle
    a.Leapfrog()
    a.Record()
    a.Update()
    n = int(round(numYears * 36.525, -1))
    for t in range(n):
        if 100*(t)/n % 20 == 0:
            print(100*(t)/n,"%", end=" ",flush=True)
        a.Record()
        a.doTimestep()
    print("100%", end=" ", flush=True)
    a.Record()
    a.File.close()

def solarSystemTests():
    m_other   = 0.33011 * 10**24 # mercury
    sem_other = 57.91   * 10**9
    # m_other   = 4.8675 * 10**24 # venus
    # sem_other = 108.21 * 10**9
    # m_other   = 0.64171 * 10**24 # mars
    # sem_other = 227.92  * 10**9
    # m_other   = 1898.19 * 10**24 # Jup
    # sem_other = 778.57  * 10**9
    # m_other   = 568.34  * 10**24 # Saturn
    # sem_other = 1433.53 * 10**9
    # m_other   = 86.813  * 10**24 # Uranus
    # sem_other = 2872.46 * 10**9
    # m_other   = 102.413 * 10**24 # Neptune
    # sem_other = 4495.06 * 10**9
    a = System(2 * 24 * 60 * 60, file="SOLARSYSTEST.bin") # 2 days
    v_other = vel(sem_other, sem_other)
    a.AddParticle(M, 0.0, 0.0, 0.0, 0.0, 0, 0.0) # Sun
    a.AddParticle(m_other, sem_other, 0.0, 0.0, 0.0, +v_other, 0.0) # Other
    v_earth = vel(sem_ear, sem_ear)
    a.AddParticle(m_earth, -sem_ear, 0, 0, 0, -v_earth, 0) # Earth
    a.Leapfrog()
    a.Record()
    a.Update()
    dur = 100 * 365.25 / 2
    for n in range(int(dur)):
        a.Record()
        a.doTimestep()
        if n % 1000 == 0:
            print(n, end=" ", flush=True)
    a.Record()
    a.File.close()
    setAnimate(widthheight_=scale[0]*1.2, scale_=scale, 
                timescale_=timescale, tracelength_=tracelength)
    animateFile("SOLARSYSTEST.bin", p=3, frameskip=frameskip, repeat=False, ax=(0,1))
    graphR("SOLARSYSTEST.bin")

def randomP():
    x,y,z,vx,vy,vz = genRandomPosVel()
    return np.array([[x,y,z], [vx,vy,vz], [0.0,0.0,0.0]])

def testEnv():
    te = System(1*24*60*60, file="TEST.bin")

    v = vel(aph_ear, aph_ear)
    v2 = vel(aph_ear, 0.75*aph_ear)

    te.AddParticle(M, 0,0,0, 0,-m_earth/2 * v / M + m_earth/2 * v2 / M,0)
    
    te.AddParticle(m_earth/2, aph_ear,0,0, 0,v,0)
    te.AddParticle(m_earth/2, -aph_ear,0,0, 0,-v2,0)
    # te.AddParticle(M, -20*eps,0,0, 0.002*eps,0,0)
    # te.AddParticle(M, -5*eps,0,0, 0,0,0)
    for x in range(365*15):
        te.Record() 
        te.doTimestep()
        #te.doTimestep(0.01*60*60)
    te.Record()
    te.File.close()

    p=3

    setAnimate(widthheight_=scale[0]*1.2, scale_=scale, 
                timescale_=timescale, tracelength_=5)
    animateFile("TEST.bin", p=p, frameskip=10)
    graphEnergies("TEST.bin", change=False, p=p)
    graphEnergies("TEST.bin", change=True, p=p)

if __name__=="__main__":
    from time import perf_counter
    from random import seed

    # seed(123456789)
    # t0 = perf_counter()
    # solarCollapseNoJup(p-2, f"SolColNoJupL10000.bin", 10000)
    # t1 = perf_counter()
    # q = t1 - t0
    # if q > 3600:
    #     print(f"No Jupiter: {q//3600}h {(q%3600)//60}m {q%60}s")
    # elif q > 60:
    #     print(f"No Jupiter: {q//60}m {q%60}s") 
    # else:
    #     print(f"No Jupiter: {q}s")

    # seed(123456789)
    # t0 = perf_counter()
    # solarCollapseJup(p-2, f"SolColJupL.bin", 10000)
    # t1 = perf_counter()
    # q = t1 - t0
    # if q > 3600:
    #     print(f"Jupiter: {q//3600}h {(q%3600)//60}m {q%60}s")
    # elif q > 60:
    #     print(f"Jupiter: {q//60}m {q%60}s")
    # else:
    #     print(f"Jupiter: {q}s")
    # rms radius and median radius -> half mass radius for glob
    # for j in range(5):
    #     t0 = perf_counter()
    #     solarCollapseJup(p-2, f"SolColJup2{j}.bin")
    #     t1 = perf_counter()
    #     q = t1 - t0
    #     if q > 3600:
    #         print(f"Jupiter {j}: {q//3600}h {(q%3600)//60}m {q%60}s")
    #     elif q > 60:
    #         print(f"Jupiter {j}: {q//60}m {q%60}s")
    #     else:
    #         print(f"Jupiter {j}: {q}s")

    # for j in range(3, 5):
    #     t0 = perf_counter()
    #     solarCollapseNoJup(p-2, f"SolColNoJup2{j}.bin")
    #     t1 = perf_counter()
    #     q = t1 - t0
    #     if q > 3600:
    #         print(f"No Jupiter {j}: {q//3600}h {(q%3600)//60}m {q%60}s")
    #     elif q > 60:
    #         print(f"No Jupiter {j}: {q//60}m {q%60}s")
    #     else:
    #         print(f"No Jupiter {j}: {q}s")

    # testEnv()
    # DEFAULTFILE = "TEST.bin"
    # p=5

    # solarSystemTests()

    setAnimate(widthheight_=scale[0]*1.2, scale_=scale, 
                timescale_=timescale, tracelength_=2, 
                shiftToFirstParticle_=False)
    # # animateFile("SolColJupL10000.bin", p=102, frameskip=1000, repeat=False, ax=(0,1))

    graphEnergies("SolColNoJupL10000.bin", False, p=101, step=1000)
    graphEnergies("SolColNoJupL10000.bin", True, p=101, step=1000)
    # animateFile("SolColNoJupL10000.bin", p=101, frameskip=10000, repeat=False, ax=(0,1))
    # animateFile(DEFAULTFILE, p=p, frameskip=1, repeat=False, ax=(0,2))
