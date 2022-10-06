from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

from Shapes import *

from random import randint, random
from math import log, sin, radians
from itertools import combinations
window = 0
width, height = 750, 750
# mouseX, mouseY = 0, 0 # unused?

G = 1          # AU
#G = 6.67 * 10**-11 # m^3 kg^-1 s^-2
stepsPerFrame = 0 # number of timesteps before each frame is displayed
timestep = 0.01    # s
t = 0             # s
arrows = True # velocity and acceleration arrows
lengthScale = 1

f = open("Energy.txt", "w") # text file used to record the energies to graph later

class Planet(Circle):
    def __init__(self, name, rad, mass, colour, x, y):
        super().__init__(x, y, rad, 20, colour)
        self.name = name
        self.mass = mass
        self.acc = Vector(0, 0)
        self.vel = Vector(0, 0)
        self.pos = Vector(x, y)
        
        self.path = []
        self.drawpath = True
        
        self.c = 99
    
    @classmethod
    def fromorbit(cls, name, rad, ownmass, parent, colour, semimajoraxis, aphelion, angle):
        pos = Vector.frompolar(aphelion, angle)
        pl = cls(name, rad, ownmass, colour, parent.pos.x + pos.x, parent.pos.y + pos.y)
        
        if semimajoraxis == 0:
            vel = Vector(0, 0)
        else:
            # the velocity of the orbit is at a right angle to the position vector
            vel = pos.normalise()
            vel = vel.rotateby(90) # degrees
            # this equation only works when the parent mass is much larger than child
            #vel *= (G * parent.mass / semimajoraxis)**0.5
            # this equation only works for all non-extreme orbits
            vel *= (G * parent.mass * (2 / aphelion - 1 / semimajoraxis)) ** 0.5
        pl.vel = vel + parent.vel
        return pl
        
    def interact(self, objects):
        for _ in range(stepsPerFrame):
            self.acc = Vector(0, 0)
            for x in objects:
                if self != x:
                    d = self.dist(x)
                    #self.acc -= d.normalise() * G * x.mass / ((abs(d)**2+3/2*ep**2))
                    self.acc -= d * G * x.mass / (abs(d)**2 + self.rad**2)**(3/2)
            self.vel += self.acc * timestep
            self.c += abs(self.vel)**2
            if self.c >= 100:
                self.c = 0
                self.path.append([self.pos.x * lengthScale, self.pos.y * lengthScale])
                if len(self.path) > 100:
                    self.path.pop(0)
            self.move(self.vel.x * timestep, self.vel.y * timestep)
            
    def draw(self):
        self.pos *= lengthScale
        super().draw()
        if self.drawpath:
            glLineWidth(1)
            glBegin(GL_LINES)
            glColor3fv(self.fillcolour)
            # each glVertex2fv _pair_ creates a single line
            for i in range(len(self.path)-1):
                glVertex2fv(self.path[i])
                glVertex2fv(self.path[i+1])
            if len(self.path) > 0:
                glVertex2fv(self.path[-1])
                glVertex2fv((self.pos.x, self.pos.y))
            glEnd()
        
        if arrows:
            vdir = self.vel * 10 * lengthScale
            glLineWidth(1)
            glColor3fv((255, 0, 255))
            glBegin(GL_LINES)
            glVertex2fv((self.pos.x, self.pos.y))
            glVertex2fv((self.pos.x + vdir.x, self.pos.y + vdir.y))          
            glEnd()
            
            glColor3fv((0, 255, 255))
            adir = self.acc * 200 * lengthScale
            glBegin(GL_LINES)
            glVertex2fv((self.pos.x, self.pos.y))
            glVertex2fv((self.pos.x + adir.x, self.pos.y + adir.y))
            glEnd()
        self.pos /= lengthScale

    def __str__(self):
        return self.name
    
    def __repr__(self):
        return f"Planet({self.name})"
    
    def kineticEnergy(self):
        # kg (ms^-1)^2 -> kg m^2 s^-2
        return 0.5 * self.mass * abs(self.vel) ** 2
    
    def gravitationalPotentialEnergy(self, obj):
        d = self.dist(obj)
        # m^3 kg^-1 s^-2 kg^2 m^-2 -> kg m s^-2
        # gpe = -G * x.mass * self.mass / (abs(d)**2)
                
        # m^3 kg^-1 s^-2 kg^2 m^-1 -> kg m^2 s^-2 
        #gpe = -G * obj.mass * self.mass / abs(d)
        gpe = -G *obj.mass * self.mass / (abs(d)**2 + self.rad**2)**0.5
        return gpe

def RealValueSystem():
    '''scaled to have the Sun be 100 units of mass, 0.07 units of radius'''
    Sun = {"Name": "Sun",
        "Mass": 1.989 * 10**30, # kg
        "Radius": 6.960 * 10**5, # km
        "Colour": (255, 255, 0),
        "Semimajor": 0, # km, average distance in orbit
        "Aphelion": 0 # km, highest point in orbit
        } # Yellow

    Mercury = {"Name": "Mercury",
            "Mass": 0.33011 * 10**24, 
            "Radius": 2439.7, 
            "Colour": (102, 0, 0),
            "Semimajor": 57.91 * 10**6,
            "Aphelion": 69.82 * 10**6
            } # Grey

    Venus = {"Name": "Venus",
            "Mass": 4.8675 * 10**24,
            "Radius": 6051.8,
            "Colour": (255, 163, 51),
            "Semimajor": 108.21 * 10**6,
            "Aphelion": 108.94 * 10**6
            } # ?Weak? Yellow
    
    Earth = {"Name": "Earth",
            "Mass": 5.9723 * 10**24,
            "Radius": 6371.000,
            "Colour": (0, 0, 255),
            "Semimajor": 149.60 * 10**6,
            "Aphelion": 152.10 * 10**6
            } # Blue

    Mars = {"Name": "Mars",
            "Mass": 0.64171 * 10**24,
            "Radius": 3389.5,
            "Colour": (255, 0, 0),
            "Semimajor": 227.92 * 10**6,
            "Aphelion": 249.23 * 10**6
            } # Red

    Jupiter = {"Name": "Jupiter",
            "Mass": 1898.19 * 10**24,
            "Radius": 69911,
            "Colour": (204, 204, 0),
            "Semimajor": 778.57 * 10**6,
            "Aphelion": 816.62 * 10**6
            } # Orange

    Saturn = {"Name": "Saturn",
            "Mass": 568.34 * 10**24,
            "Radius": 58232,
            "Colour": (255, 255, 102),
            "Semimajor": 1433.53 * 10**6,
            "Aphelion": 1514.50 * 10**6
            } # Yellow

    Uranus = {"Name": "Uranus",
            "Mass": 86.813 * 10**24,
            "Radius": 25362,
            "Colour": (0, 255, 255),
            "Semimajor": 2872.46 * 10**6,
            "Aphelion": 3003.62 * 10**6
            } # light blue

    Neptune = {"Name": "Neptune",
            "Mass": 102.413 * 10**24,
            "Radius": 24622,
            "Colour": (0, 0, 255),
            "Semimajor": 4495.06 * 10**6,
            "Aphelion": 4545.67 * 10**6
            } # blue 

    Planets = [Mercury, Venus, Earth, Mars, 
            Jupiter, Saturn, Uranus, Neptune]
    
    #massScale = 100 / Sun["Mass"]
    global lengthScale
    lengthScale = 700 / Neptune["Aphelion"]

    system = [Planet(Sun["Name"], 2, Sun["Mass"], Sun["Colour"], width/2, height/2)]

    for obj in Planets:
        p = Planet.fromorbit(
            obj["Name"],                       # Name
            2, # obj["Radius"]# * lengthScale,       # Radius
            obj["Mass"], #* massScale,           # Mass
            system[0],                         # Parent, Sun
            obj["Colour"],                     # Colour
            obj["Semimajor"],# * lengthScale,    # Semimajor axis
            obj["Aphelion"],# * lengthScale,     # Aphelion
            randint(0, 360)                    # (random) angle
        )
        #a = obj["Semimajor"], obj["Aphelion"]
        system.append(p)
    
    EarthMoon = {"Name": "EarthMoon",
            "Mass": 7.346 * 10**22,
            "Radius": 1737.4,
            "Colour": (128, 128, 128),
            "Semimajor": 384400,
            "Aphelion": 405500,
            "Parent": system[3]
            } # Grey
    
    # Moons = [EarthMoon]
    
    # for moon in Moons:
    #     p = Planet.fromorbit(
    #         moon["Name"],
    #         moon["Radius"] * lengthScale,
    #         moon["Mass"] * massScale,
    #         moon["Parent"],
    #         moon["Colour"],
    #         moon["Semimajor"] * lengthScale,
    #         moon["Aphelion"] * lengthScale,
    #         randint(0, 360)
    #     )
    #     system.append(p)
    
    return system

def randomSystem(n):
    return [Planet(str(q), 5, 5, (random(),random(),random()), randint(0, width), randint(0,height)) for q in range(n)]

#system = RealValueSystem()

# system = [Planet("SUN", 5, 100, (100,100,100), width/2 , height/2)]
# system.append(Planet.fromorbit(name="Alpha", rad=2.5, ownmass=0.1, 
#                                 parent=system[0], colour=(200, 200, 0), 
#                                 semimajoraxis=100, aphelion=150, angle=0))


# bary = Planet("BARYCENTER", 1, 100, (0,0,0), width/2, height/2)
# system = [Planet.fromorbit("STAR A", 5, 100, bary, (1, 0, 0), 50, 75, 0),
#           Planet.fromorbit("STAR B", 5, 100, bary, (0, 0, 1), 50, 75, 180)]

system = randomSystem(10)

def refresh2d(width, height):
    #glViewport(int(width/2-system[0].pos.x), int(height/2-system[0].pos.y), width, height)
    glViewport(0, 0, width, height)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    glOrtho(0, width, 0, height, 0, 1)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    
def draw():
    global t
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glLoadIdentity()
    refresh2d(width, height)
    
    frame = system[:]
    
    for x in system:
        x.interact(frame)
        x.draw()
        
    del frame
    t += timestep * stepsPerFrame
    overlay()    
    
    glutSwapBuffers()    

def handleKeyinput(key: bin, x: int, y: int):
    global arrows
    global stepsPerFrame
    global timestep
    if key == b"1":
        # toggle vel and acc arrows
        arrows = not arrows
        print(f"Arrows: {arrows}")
    elif key == b"2":
        # increase steps per frame
        if stepsPerFrame <= 0:
            stepsPerFrame = 1
        else:
            stepsPerFrame *= 10
        print(f"stepsPerFrame: {stepsPerFrame}")
    elif key == b"3":
        # decrease steps per frame
        if stepsPerFrame <= 1:
            stepsPerFrame = 0
        else:
            stepsPerFrame /= 10
        stepsPerFrame = int(stepsPerFrame)
        print(f"stepsPerFrame: {stepsPerFrame}")
    elif key == b"4":
        # increase timestep
        timestep *= 10
        print(f"timestep: {timestep}")
    elif key == b"5":
        # decrease timestep
        timestep /= 10
        print(f"timestep: {timestep}")
    elif key == b"6":
        timestep *= -1
        print(f"timestep: {timestep}")
    else:
        print(f"Key {key} has no mapping")
    
def overlay():
    glMatrixMode(GL_PROJECTION)
    glPushMatrix()
    glLoadIdentity()
    
    glOrtho(0, width, 0, height, -1, 1)
    
    glMatrixMode(GL_MODELVIEW)
    glPushMatrix()
    glLoadIdentity()
    
    glColor3fv((255, 0, 0))
    
    glDisable(GL_DEPTH_TEST)
    #glDisable(GL_LIGHTING)
    
    ke = getKineticEnergy()
    gpe = getGPE()
    
    f.write(str(ke) + "," + str(gpe) + "," + str(gpe + ke) + "," + str(stepsPerFrame * timestep) + "\n")
    
    
    menu = f"""Arrows: {arrows}
Toggle with 1
StepsPerFrame: {stepsPerFrame}
Inc with 2; Dec with 3
TimeStep: {timestep}
Inc with 4; Dec with 5; Invert with 6
Kinetic Energy: {ke}
Gravitational Potential Energy: {gpe}
Total Energy: {ke + gpe}
G: {G}"""
    
    menu = menu.split("\n")
    for string, x in zip(menu, range(10 + 15 * len(menu), 10, -15)):
        glRasterPos2i(10, x)
        for c in string:
            glutBitmapCharacter(GLUT_BITMAP_8_BY_13, ord(c))
            
    glEnable(GL_DEPTH_TEST)
    #glEnable(GL_LIGHTING)

    glMatrixMode(GL_MODELVIEW)
    glPopMatrix()
    glMatrixMode(GL_PROJECTION)
    glPopMatrix()

def getKineticEnergy():
    return sum(map(lambda x: x.kineticEnergy(), system))

def getGPE():
    return sum(map(lambda x: x[0].gravitationalPotentialEnergy(x[1]), combinations(system, 2)))

if __name__ == "__main__":
    glutInit()
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH)
    glutInitWindowSize(width, height)
    glutInitWindowPosition(25, 25)
    glutCreateWindow(b"Title")
    glutDisplayFunc(draw)
    glutIdleFunc(draw)
    glutKeyboardFunc(handleKeyinput)
    # glutOverlayDisplayFunc(overlay)
    # glutShowOverlay()
    glutMainLoop()
    
    f.close()
