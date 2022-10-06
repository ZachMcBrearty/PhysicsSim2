from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

from Shapes import Circle, Vector

from math import floor, atan, cos, sin, radians
import colorsys

chargescale = 1.0

q = 1.6022*10**-19 # C
pi=3.1415
ep0=8.85*10**-12
k = 1/(4*pi*ep0)

class Charge:
    def __init__(self, charge, x, y):
        self.charge = charge
        self.x = x
        self.y = y
        self.pos = Vector(x, y)

    def __str__(self):
        return f"Charge({self.charge}, {self.x}, {self.y})"
    __repr__=__str__
    def fieldAtPoint(self, x, y):
        if x == self.x and y == self.y:
            return 0
        else:
            return k * self.charge / ((x - self.x)**2 + (y - self.y)**2)
    
    def fieldAtPointVect(self, v):
        dir = v - self.pos
        if abs(dir)==0:
            return Vector(0, 0)
        # k * q / r**2 * unit vector from charge to point
        return k * self.charge / abs(dir)**2 * dir.normalise()
    
class ScalarElectricField:
    def __init__(self, x, y, width, height, density, charges=[]):
        '''a flat plane with the z component shown as colour'''
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        self.density = density
        
        w = floor(width/density)
        h = floor(height/density)
        self.field = [[0 for a in range(w)] 
                    for b in range(h)]
        self.charges = charges
        if charges != []:
            self.calculateField()
    
    def addCharge(self, charge):
        assert(isinstance(charge, Charge))
        self.charges.append(charge)
        self.calculateField()
     
    def addCharges(self, charges):
        #assert(isinstance(charge, Charge))
        self.charges.extend(charges)
        self.calculateField()
        
    def calculateField(self):
        for charge in self.charges:
            for x in range(len(self.field)):
                for y in range(len(self.field[x])):
                    self.field[x][y] += charge.fieldAtPoint(x, y)
        
    def draw(self):
        glPointSize(1)
        glBegin(GL_POINTS)
        
        scale = atan(255/255)
        for i, v in enumerate(self.field):
            for j, w in enumerate(v):
                if round(w*chargescale, 1) == 0:
                    glColor3fv((1,1,1))
                    glVertex2fv((self.x + i * self.density, self.y + j * self.density))
                else:
                    h = ((atan(w * chargescale) / scale) + 1)/2
                    c = colorsys.hls_to_rgb(h, 0.5, 0.5)
                    glColor3fv(c)
                    glVertex2fv((self.x + i * self.density, self.y + j * self.density))

        glEnd()
        glPointSize(4)
        glBegin(GL_POINTS)
        glColor3fv((1,1,1))
        for charge in self.charges:
            glVertex2fv((self.x + charge.x * self.density, self.y + charge.y * self.density))

        glEnd()

class VectorElectricField:
    def __init__(self, x, y, width, height, density, charges=[]):
        '''a flat plane with the z component shown as colour'''
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        self.density = density
        
        w = floor(width/density)
        h = floor(height/density)
        self.field = [[Vector(0,0) for a in range(w)] 
                    for b in range(h)]
        self.charges = charges
        if charges != []:
            self.calculateField()
    
    def addCharge(self, charge):
        assert(isinstance(charge, Charge))
        self.charges.append(charge)
        self.calculateField()
     
    def addCharges(self, charges):
        #assert(isinstance(charge, Charge))
        self.charges.extend(charges)
        self.calculateField()
        
    def calculateField(self):
        for x in range(len(self.field)):
            for y in range(len(self.field[x])):
                e = Vector(0,0)
                for charge in self.charges:
                    e += charge.fieldAtPointVect(Vector(x, y))
                self.field[x][y] += e
        
    def draw(self):
        glPointSize(1)
        glBegin(GL_POINTS)
        
        scale = atan(255/255)
        for i, v in enumerate(self.field):
            for j, w in enumerate(v):
                if round(abs(w*chargescale), 1) == 0:
                    glColor3fv((1,1,1))
                    glVertex2fv((self.x + i * self.density, self.y + j * self.density))
                else:
                    h = ((atan(abs(w*chargescale)) / scale) + 1)/2
                    c = colorsys.hls_to_rgb(h, 0.5, 0.5)
                    glColor3fv(c)
                    glVertex2fv((self.x + i * self.density, self.y + j * self.density))

        glEnd()
        glPointSize(4)
        glBegin(GL_POINTS)
        glColor3fv((1,1,1))
        for charge in self.charges:
            glVertex2fv((self.x + charge.x * self.density, self.y + charge.y * self.density))

        glEnd()

width, height = 1000, 1000
dens = 5
v = VectorElectricField(50, 50, width-100, height-100, dens)
#c = Circle(500-5/7*18*dens, 500-5/14*18*dens, 5*3**0.5/7 * 18*dens, 1, strokecolour=(1,1,1), strokewidth=0.1)
### create charges
# q1 = Charge(-2.5*10**-6, 90, 90)
# q2 = Charge(6*10**-6, 108, 99)
# v.addCharges([q1, q2])
# q1 = Charge(2, 90, 90)
# q2 = Charge(-1, 108, 90)
# q3 = Charge(-1, 90-18, 90)
# v.addCharges([q1, q2, q3])
q1 = Charge(1, 90-18, 90)
q2 = Charge(-4, 108, 90)
v.addCharges([q1, q2])
# charges = []
# R = 18
# for x in range(0, 6):
#     if x != 3:
#         q1 = Charge(1, 90+R*cos(radians(60*x+90)), 90+R*sin(radians(60*x+90)))
#         charges.append(q1)
# print(charges)
# v.addCharges(charges)

def refresh2d(width, height):
    glViewport(0, 0, width, height)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    glOrtho(0, width, 0, height, 0, 1)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    
def draw():
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glLoadIdentity()
    refresh2d(width, height)
    
    #v.interact()
    v.draw()
    #Circle(-5/7, -5/14, 5*3**0.5/7, 1)
    #c.draw()
    
    glutSwapBuffers()

def d():pass

def handleKeyinput(key: int, x: int, y: int):
    global chargescale
    if key == b"1":
        # inc charge
        chargescale *= 10
        print(f"chargescale: {chargescale}")
    elif key == b"2":
        # dec charge
        chargescale /= 10
        print(f"chargescale: {chargescale}")
    else:
        print(f"Key {key} has no mapping")
    draw()
    
if __name__ == "__main__":
# mypy: ignore-errors
    glutInit()
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH)
    glutInitWindowSize(width, height)
    glutInitWindowPosition(0, 0)
    glutCreateWindow(b"Title")
    glutDisplayFunc(d)
    glutIdleFunc(d)
    glutKeyboardFunc(handleKeyinput)
    draw()
    draw()
    glutMainLoop()