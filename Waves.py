from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

from Shapes import *
from math import floor, atan

import colorsys

timestep = 0.1
stepsPerFrame = 1
smooth = True

class Vector3DField:
    def __init__(self, x, y, width, height, density, wavespeed, bound=True):
        '''a flat plane with the z component shown as colour'''
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        self.density = density
        self.wavespeed = wavespeed
        
        self.boundary = bound
        w = floor(width/density)
        h = floor(height/density)
        self.pos = [[0 for a in range(w)] 
                    for b in range(h)]
        self.vel = [[0 for a in range(w)] 
                    for b in range(h)]
        self.acc = [[0 for a in range(w)] 
                    for b in range(h)]
    
    def lapacian(self, i, j):
        if self.boundary:
            if (i == 0 or i == len(self.pos)-1 or j == 0 or j == len(self.pos[0])-1):
                return 0
            else:
                s = 0
                c = 0
                for x in range(-1, 2):
                    for y in range(-1, 2):
                        s += self.pos[i + x][j + y]
                        c += 1
                psi_av = s / c
                psi = self.pos[i][j]
                
                return  (6 / self.density**2) * (psi_av - psi)
        else:
            s = 0
            c = 0
            for x in range(-1, 2):
                if i + x >= len(self.pos) or i + x < 0:
                    continue
                for y in range(-1, 2):
                    if j + y >= len(self.pos[i+x]) or j + y < 0:
                        continue
                    s += self.pos[i + x][j + y]
                    c += 1
            psi_av = s / c
            psi = self.pos[i][j]
        
            return  (6 / self.density**2) * (psi_av - psi)
    
    def interact(self):
         # acc = c**2 * lap(psi)
        for _ in range(stepsPerFrame):
            for i in range(len(self.pos)):
                for j in range(len(self.pos[i])):
                    self.acc[i][j] = (self.wavespeed**2) * self.lapacian(i, j)
        
            for i in range(len(self.pos)):
                for j in range(len(self.pos[i])):
                    self.vel[i][j] += self.acc[i][j] * timestep
                    self.pos[i][j] += self.vel[i][j] * timestep
        
    def draw(self):
        glBegin(GL_POINTS)
        
        scale = atan(255/255)
        
        for i, v in enumerate(self.pos):
            for j, w in enumerate(v):
                if smooth:
                    h = ((atan(w / 255) / scale) + 1)/2
                else:
                    h = round(((atan(w / 255) / scale) + 1)/2, 1)
                c = colorsys.hls_to_rgb(h, 0.5, 0.5)
                glColor3fv(c)
                glVertex2fv((self.x + i * self.density, self.y + j * self.density))

        glEnd()
        
class VectorField:
    def __init__(self, x, y, width, height, a, c):
        self.wavespeed = c
        self.a = a
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        self.pos = [Vector(b,0) for b in range(0, width+1, a)]
        self.vel = [Vector(0,0) for b in range(0, width+1, a)]
        self.acc = [Vector(0,0) for b in range(0, width+1, a)]
        
    def lapacian(self, i):
        # lap(psi) = 6/a**2 (psi_av - psi)
        if i == 0:
            # psi_av = (self.pos[0].y + self.pos[1].y) / 2
            psi_av = self.pos[0].y
        elif i == len(self.pos) - 1:
            # psi_av = (self.pos[-1].y + self.pos[-2].y) / 2
            psi_av = self.pos[-1].y
        else:
            psi_av = (self.pos[i-1].y + self.pos[i].y + self.pos[i+1].y) / 3

        return (6 / self.a**2) * (psi_av - self.pos[i].y)
    
    def interact(self):
        '''attempt to "flatten" the field with the lapacian'''
        # acc = c**2 * lap(psi)
        for _ in range(stepsPerFrame):
            for i in range(len(self.pos)):
                self.acc[i].y = (self.wavespeed**2) * self.lapacian(i)
        
            for i in range(len(self.pos)):
                self.vel[i] += self.acc[i] * timestep
                self.pos[i] += self.vel[i] * timestep
        
    def draw(self):
        glBegin(GL_POINTS)
        glColor3f(1, 1, 1)
        
        for v in self.pos:
            glVertex2fv((self.x + v.x, self.y + v.y))

        glEnd()
        glBegin(GL_LINES)
        glColor3f(1,2,1)
        
        for i in range(len(self.pos)-1):
            glVertex2fv((self.x + self.pos[i].x, self.y + self.pos[i].y))
            glVertex2fv((self.x + self.pos[i+1].x, self.y + self.pos[i+1].y))
        glEnd()

def perturb(VecField, n, size):
    step = n * 180 / (VecField.width / VecField.a)
    for i in range(len(VecField.pos)):
        VecField.pos[i].y += size * sin(radians(i * step))

def drive(VecField):
    mid = len(VecField.pos) // 2
    for x in range(0, 91, 1):
        VecField.pos[mid - 45+x].y = 250 * sin(2*radians(x))

def perturb3d(Vec3DField, n, m, size):
    stepX = n * 180 / len(Vec3DField.pos)
    stepY = m * 180 / len(Vec3DField.pos[0])
    for i in range(len(Vec3DField.pos)):
        for j in range(len(Vec3DField.pos[i])):
            
            Vec3DField.pos[i][j] += size * sin(radians((i+0.5) * stepX)) * sin(radians((j+0.5) * stepY))
            
def drive3d(Vec3DField, width, height):
    midX = len(Vec3DField.pos) // 2
    midY = len(Vec3DField.pos[0]) // 2
    for x in range(0, width+1, 1):
        for y in range(0, height+1):
            Vec3DField.pos[midX - round(width/2)+x][midY-round(height/2)+y] += 250 * sin(180/width*radians(x)) * sin(radians(180/height*y))
    
width, height = 500, 500
# v = VectorField(50, height/2, width-100, height-100, 1, 5)
# #drive(v)
# for x in range(1, 2, 2):
#     perturb(v, x, 50)
dens = 10
v = Vector3DField(50, 50, width-100, height-100, dens, 10)
perturb3d(v, 1, 1, 512)
#drive3d(v, int((width-100)/(2*dens)), int((height-100)/(2*dens)))
#drive3d(v, 10, 90)

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
    
    v.interact()
    v.draw()
    
    overlay()
    
    glutSwapBuffers()    

def handleKeyinput(key: bin, x: int, y: int):
    global smooth
    global stepsPerFrame
    if key == b"1":
        # increase steps per frame
        if stepsPerFrame <= 0:
            stepsPerFrame = 1
        else:
            stepsPerFrame *= 10
        print(f"stepsPerFrame: {stepsPerFrame}")
    elif key == b"2":
        # decrease steps per frame
        if stepsPerFrame <= 1:
            stepsPerFrame = 0
        else:
            stepsPerFrame /= 10
        stepsPerFrame = int(stepsPerFrame)
        print(f"stepsPerFrame: {stepsPerFrame}")
    elif key == b"3":
        smooth = not smooth
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
    
    menu = f"""StepsPerFrame: {stepsPerFrame}; Inc with 1; Dec with 2
Smooth Colour: {smooth}; Toggle with 3"""
    
    menu = menu.split("\n")
    for string, x in zip(menu, range(10 + 15 * (len(menu)-1), 9, -15)):
        glRasterPos2i(10, x)
        for c in string:
            glutBitmapCharacter(GLUT_BITMAP_8_BY_13, ord(c))
            
    glEnable(GL_DEPTH_TEST)

    glMatrixMode(GL_MODELVIEW)
    glPopMatrix()
    glMatrixMode(GL_PROJECTION)
    glPopMatrix()

if __name__ == "__main__":
    print(bool(glutInit))
    glutInit()
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH)
    glutInitWindowSize(width, height)
    glutInitWindowPosition(0, 0)
    glutCreateWindow(b"Title")
    glutDisplayFunc(draw)
    glutIdleFunc(draw)
    glutKeyboardFunc(handleKeyinput)
    glutMainLoop()