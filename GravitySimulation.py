from itertools import combinations
from more_itertools import flatten
import numpy as np
from io import TextIOWrapper as _TextIOWrapper

G = 1
scale = 100
timestepB = 100

eps = 0.1

DEFAULTFILE = "GravitySim.bin"

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
        self.minDist = scale

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
        
    
    def Update(self) -> None:
        """Update the positions and velocities of the particles in the system based
        on the current time step"""
        self.CurTimeStep = self.MaxTimeStep * (np.log(self.minDist/timestepB + 1) / np.log(scale/timestepB + 1) + .001) 
        updateMatrix = np.array([[1, self.CurTimeStep, 0],[0, 1, self.CurTimeStep],[0, 0, 0]])
        self.Particles = updateMatrix @ self.Particles
        self.minDist = scale
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

if __name__=="__main__":
    # ener = []
    # a = System(0.01)
    # a.AddParticle(2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    # a.AddParticle(2.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    # #print(a.Particles[:,:,0])
    # for _ in range(1000):
    #     a.doTimestep(2)
    #     a.Record()
    #     #print(a.Particles[:,:,0])
    #     k, g = a.kineticEnergy(), a.gravitationalEnergy()
    #     ener.append(k+g)
    # ener = np.array(ener)
    # print(list(100*abs(1 - ener / ener[0])))
    # a.File.close()
    
    # Example Reading of the data file
    binfile = open(DEFAULTFILE, "rb")
    b = np.fromfile(binfile)
    n = 3 + 3 * 2 # 3 -> time, kinetic, gravity, 3 -> 3d, 2 = number of particles
    #print("TIME : ", b[0::n])
    #print("KINET: ", b[1::n])
    #print("GRAVI: ", b[2::n])
    #print("TOTAL: ", b[1::n]+b[2::n])
    t0 = (b[1::n]+b[2::n])[0]
    te = (b[1::n]+b[2::n])[-1]
    print(100*abs(1-te/t0))
    print(np.average(b[3::n]))
    print(np.std(b[3::n]))
    print(np.average(b[6::n]))
    print(np.std(b[6::n]))
    #print("PART1: ", b[3::n], b[4::n], b[5::n])
    #print("PART2: ", b[6::n], b[7::n], b[8::n])
    binfile.close()