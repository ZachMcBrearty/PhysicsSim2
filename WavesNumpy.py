from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

import numpy
from numpy import floor, arctan as atan, sin, radians, round, exp
from numpy import cos, pi, ceil

from random import random

import colorsys

setup = True
setupStage = 0
_2d=False
selection = -1
paramsNum = 0
params = [0,0,0,0,0]
v = None

unbound = False
bound = True
continuous = False

timestep = 0.1
stepsPerFrame = 0
smooth = True
velocities = True
accelerations = True

bytesNums = [bytes(str(x), "utf8") for x in range(0, 10)]

width, height = 800, 800

class Vector2DField:
    def __init__(self, x, y, width, density, wavespeed):
        self.x = x
        self.y = y
        self.width = width
        self.density = density
        self.wavespeed = wavespeed

        self.w = int(floor(width/density))

        self.pos = numpy.zeros(shape=(self.w))
        self.vel = numpy.zeros(shape=(self.w))
        self.acc = numpy.zeros(shape=(self.w))

    def lapacian(self):
        lap = numpy.zeros(shape=self.pos.shape)
        for i in range(1, self.w-1):
            avg = numpy.average(self.pos[i-1:i+2])
            psi = self.pos[i]
            lap[i] = 6 * (avg - psi) / self.density**2
            lap[i] *= 1+(random()-0.5)/10
        if bound:
            pass
        elif unbound:
            lap[0] = 6 * (self.pos[1] - self.pos[0]) / self.density ** 2
            lap[-1] = 6 * (self.pos[-2] - self.pos[-1]) / self.density ** 2
        else: # continuous
            lap[0] = 6 * ((self.pos[-1] + self.pos[0] + self.pos[1])/3 - self.pos[0]) / self.density**2
            lap[-1] = 6 * ((self.pos[-2] + self.pos[-1]+ self.pos[0])/3 - self.pos[-1]) / self.density**2
        return lap

    def interact(self):
        for _ in range(stepsPerFrame):
            self.acc = self.wavespeed**2 * self.lapacian()
            
            self.vel += self.acc * timestep
            self.pos += self.vel * timestep
    
    def draw(self):
        if smooth:
            glColor3f(1, 0, 0)
            glBegin(GL_LINES)  
            for i in range(self.w-1):
                glVertex2f(self.x + i*self.density, self.y + self.pos[i])
                glVertex2f(self.x + (i+1)*self.density, self.y + self.pos[i+1])
            glEnd()
            if velocities:
                glColor3f(0, 1, 0)
                glBegin(GL_LINES)
                for i in range(self.w-1):
                    glVertex2f(self.x + i*self.density, self.y + self.vel[i])
                    glVertex2f(self.x + (i+1)*self.density, self.y + self.vel[i+1])
                glEnd()
            if accelerations:
                glColor3f(0, 0, 1)
                glBegin(GL_LINES)
                for i in range(self.w-1):
                    glVertex2f(self.x + i*self.density, self.y + self.acc[i])
                    glVertex2f(self.x + (i+1)*self.density, self.y + self.acc[i+1])
                glEnd()
        else:         
            glBegin(GL_POINTS)
            for i, y in enumerate(self.pos):
                glVertex2f(self.x + i*self.density, self.y + y)
            glEnd()

class Vector3DField:
    def __init__(self, x, y, width, height, density, wavespeed, bound=True):
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        self.density = density
        self.wavespeed = wavespeed

        self.bound = bound
        self.w = int(floor(width/density))
        self.h = int(floor(height/density))
        p = (self.w, self.h)
        self.pos = numpy.zeros(shape=p)
        self.vel = numpy.zeros(shape=p)
        self.acc = numpy.zeros(shape=p)
    
    def lapacian(self):
        lap = numpy.zeros(shape=self.pos.shape)
        if self.bound:
            for i in range(1, self.w-1):
                for j in range(1, self.h-1):
                    avg = numpy.average(self.pos[i-1:i+2,j-1:j+2])
                    psi = self.pos[i, j]
                    lap[i,j] = 6 / self.density**2 * (avg - psi)
                    lap[i,j] *= 1+(random()-0.5)/10
        else:
            pass
        return lap

    def interact(self):
        for _ in range(stepsPerFrame):
            self.acc = self.wavespeed**2 * self.lapacian()
            self.vel += self.acc * timestep
            self.pos += self.vel * timestep
    
    def draw(self):
        scale = atan(255/255)
        if smooth:
            h = (atan(self.pos/255) / scale + 1) / 2
        else:
            h = round((atan(self.pos/255) / scale + 1) / 2, 1)
                
        glBegin(GL_POINTS)
        for i in range(self.pos.shape[0]):
            for j in range(self.pos.shape[1]):
                col = colorsys.hls_to_rgb(h[i,j], 0.5, 0.5)
                glColor3fv(col)
                glVertex2f(self.x + i*self.density, self.y + j*self.density)
        glEnd()

def perturbStanding2DPos(vec2df, n, size):
    stepX = n*180/vec2df.w
    f = lambda i: size * sin(radians((i+0.5)*stepX))
    vec2df.pos += numpy.fromfunction(f, vec2df.pos.shape)
    
def perturbStanding2DVel(vec2df, n, size):
    stepX = n*180/vec2df.w
    f = lambda i: size * sin(radians((i+0.5)*stepX))
    vec2df.vel += numpy.fromfunction(f, vec2df.vel.shape)
    vec2df.vel[0] = 0
    vec2df.vel[-1] = 0

def perturbLocal2DPos(vec2df, x, wid, n, size, func):
    for i in range(0, wid):
        vec2df.pos[x-wid//2+i] += size * func(n, wid, i)
def perturbLocal2DVel(vec2df, x, wid, n, size, func):
    for i in range(0, wid):
        vec2df.vel[x-wid//2+i] += size * func(n, wid, i)

def sinLocal2DPos(vec2df, x, wid, n, size):
    for i in range(0, wid):
        vec2df.pos[x-wid//2+i] += size * sin(radians(n*180/wid*i))  
    vec2df.pos[0] = 0
    vec2df.pos[-1] = 0

def sinLocal2DVel(vec2df, x, wid, n, size):
    for i in range(0, wid):
        vec2df.vel[x-wid//2+i] += size * sin(radians(n*180/wid*i))
    vec2df.vel[0] = 0
    vec2df.vel[-1] = 0 

def gaussian2DPos(vec2df, x, wid, size):
    f = lambda p: size * exp(-((p - x) / wid)**2 / 2)
    vec2df.pos += numpy.fromfunction(f, vec2df.pos.shape)
    if unbound or continuous:
        vec2df.pos[0] = 0
        vec2df.pos[-1] = 0

def gaussian2DVel(vec2df, x, wid, size):
    f = lambda p: size * exp(-((p - x) / wid)**2 / 2)
    vec2df.vel += numpy.fromfunction(f, vec2df.vel.shape)
    if unbound or continuous:
        vec2df.vel[0] = 0
        vec2df.vel[-1] = 0

def func2DPos(vec2df, x, wid, size, func):
    f = lambda p: size * func(x, wid, p)
    vec2df.pos += numpy.fromfunction(f, vec2df.pos.shape)
    vec2df.pos[0] = 0
    vec2df.pos[-1] = 0

def func2DVel(vec2df, x, wid, size, func):
    f = lambda p: size * func(x, wid, p)
    vec2df.vel += numpy.fromfunction(f, vec2df.vel.shape)
    vec2df.vel[0] = 0
    vec2df.vel[-1] = 0

def perturbStanding3DPos(vec3df, n, m, size):
    stepX = n*180 / vec3df.w
    stepY = m*180 / vec3df.h
    f = lambda i, j: size * sin(radians((i+0.5)*stepX)) * sin(radians((j+0.5)*stepY))           
    vec3df.pos += numpy.fromfunction(f, vec3df.pos.shape)

def perturbStanding3DVel(vec3df, n, m, size):
    stepX = n*180 / vec3df.w
    stepY = m*180 / vec3df.h
    f = lambda i, j: size * sin(radians((i+0.5)*stepX)) * sin(radians((j+0.5)*stepY))           
    vec3df.vel += numpy.fromfunction(f, vec3df.vel.shape)

def perturbLocal3DPos(vec3df, x, y, wid, hei, size):
    for i in range(0, wid):
        for j in range(0, hei):
            vec3df.pos[x-wid//2+i, y-hei//2+j] += size * sin(radians(180/wid*i)) * sin(radians(180/hei*j))
    
def perturbLocal3DVel(vec3df, x, y, wid, hei, size):
    for i in range(0, wid):
        for j in range(0, hei):
            vec3df.vel[x-wid//2+i, y-hei//2+j] += size * sin(radians(180/wid*i)) * sin(radians(180/hei*j))
    

def refresh2d(width, height):
    glViewport(0, 0, width, height)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    glOrtho(0, width, 0, height, 0, 1)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()

def handleNumInput(key: bin, paramMax: int, func):
    global setup, setupStage, selection, _2d, v, params, paramsNum
    if key == b"\r":
        paramsNum += 1
        if paramsNum > paramMax-1:
            func(v, *params[:paramMax])
            # reset
            params = [0,0,0,0,0]
            paramsNum = 0
            selection = -1
            setupStage = 1
    elif key == b"-":
        params[paramsNum] *= -1
    elif key in bytesNums:
        params[paramsNum] *= 10
        params[paramsNum] += int(key)
    elif key == b"\x08":
        params[paramsNum] -= params[paramsNum] % 10
        params[paramsNum] =  int(params[paramsNum] / 10)
    else:
        print(f"key {key} not recognised")

def handleSetupKeyInput(key: bin, x: int, y: int):
    global setup, setupStage, selection, _2d, v, params, paramsNum, bound, unbound, continuous
    if setupStage == 0:
        if key == b"1":
            _2d = True
            v = Vector2DField(x=width/10, y=height/2, width=width*4/5, density=1, wavespeed=5)
            setupStage = 1
        elif key == b"2":
            _2d = False
            v = Vector3DField(x=width/10, y=height/10, width=width*4/5, 
                              height=height*4/5, density=10, wavespeed=10)
            setupStage = 1
        else:
            print(f"key {key} not recognised")
    elif setupStage == 1 and _2d:
        if key == b"0":
            # go back
            setup = False
        elif key in [bytes(str(x), "utf8") for x in range(1, 7)]:
            selection = int(key)
            setupStage = 2
        elif key == b"B" or key == b"b":
            unbound, bound, continuous = continuous, unbound, bound
        else:
            print(f"key {key} not recognised")   
    elif setupStage == 2 and _2d:
        if selection == 1:
            handleNumInput(key, 2, perturbStanding2DPos)
        elif selection == 2:
            handleNumInput(key, 2, perturbStanding2DVel)
        elif selection == 3:
            handleNumInput(key, 4, sinLocal2DPos)
        elif selection == 4:
            handleNumInput(key, 4, sinLocal2DVel)  
        elif selection == 5:
            handleNumInput(key, 3, gaussian2DPos)  
        elif selection == 6:
            handleNumInput(key, 3, gaussian2DVel)   

    elif setupStage == 1 and not _2d:
        if key == b"0":
            # go back
            setup = False
        elif key in [bytes(str(x), "utf8") for x in range(1, 5)]:
            selection = int(key)
            setupStage = 2
        else:
            print(f"key {key} not recognised")

    elif setupStage == 2 and not _2d:
        if selection == 1:
            handleNumInput(key, 3, perturbStanding3DPos)
        elif selection == 2:
            handleNumInput(key, 3, perturbStanding3DVel)
        elif selection == 3:
            handleNumInput(key, 5, perturbLocal3DPos)
        elif selection == 4:
            handleNumInput(key, 5, perturbLocal3DVel) 

def handleActiveKeyinput(key: bin, x: int, y: int):
    global smooth, velocities, accelerations, stepsPerFrame
    if key == b"0":
        #reset
        global setup, setupStage, v
        setup = True
        setupStage = 0
        v = None
        stepsPerFrame = 0
    elif key == b"1":
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
    elif key == b"4":
        velocities = not velocities
    elif key == b"5":
        accelerations = not accelerations
    else:
        print(f"Key {key} has no mapping")

def handleKeyInput(key: bin, x: int, y: int):
    global setup
    if setup:
        handleSetupKeyInput(key, x, y)
    else:
        handleActiveKeyinput(key, x, y)

def setupOverlay():
    if setupStage == 0:
        # 2d or 3d
        menu = f"""1. 2D wave system (String)
2. 3D wave system (Colour Map)"""
    elif setupStage == 1 and _2d:
        #Add waves to system
        menu = f"""0. Begin
1. Add Positional Standing Wave
2. Add Velocity Standing Wave
3. Add Positional Local Sine Wave
4. Add Velocity Local Sine Wave
5. Add Positional Gaussian Wave
6. Add Velocity Gaussian Wave
B. Cycle: {"*Unbound*" if unbound else "Unbound"} {"*Bound*" if bound else "Bound"} {"*Continuous*" if continuous else "Continuous"}"""
    elif setupStage == 2 and _2d:
        if selection == 1 or selection == 2:
            menu = f"""Perturb Standing {"Pos" if selection == 1 else "Vel"}
n : {params[0]}{"_" if paramsNum == 0 else ""}
Size : {params[1]}{"_" if paramsNum == 1 else ""}
Enter to Continue"""
        elif selection == 3 or selection == 4:
            menu = f"""Perturb Local Sine {"Pos" if selection == 3 else "Vel"}
x pos (0 to {v.w}) : {params[0]}{"_" if paramsNum == 0 else ""}
Width : {params[1]}{"_" if paramsNum == 1 else ""}
n : {params[2]}{"_" if paramsNum == 2 else ""}
Size : {params[3]}{"_" if paramsNum == 3 else ""}
Enter to Continue"""
        elif selection == 5 or selection == 6:
            menu = f"""Perturb Gaussian {"Pos" if selection == 5 else "Vel"}
x pos (0 to {v.w}) : {params[0]}{"_" if paramsNum == 0 else ""}
Width : {params[1]}{"_" if paramsNum == 1 else ""}
Size : {params[2]}{"_" if paramsNum == 2 else ""}
Enter to Continue"""
    elif setupStage == 1 and not _2d:
        menu = f"""0. Go Back
1. Add Standing Positional Wave
2. Add Standing Velocity Wave
3. Add Local Positional Sine Wave
4. Add Local Velocity Sine Wave"""
    elif setupStage == 2 and not _2d:
        if selection == 1 or selection == 2:
            menu = f"""Perturb Standing {"Pos" if selection == 1 else "Vel"}
n : {params[0]}{"_" if paramsNum == 0 else ""}
m : {params[1]}{"_" if paramsNum == 1 else ""}
Size : {params[2]}{"_" if paramsNum == 2 else ""}
Enter to Continue"""
        elif selection == 3 or selection == 4:
            menu = f"""Perturb Standing {"Pos" if selection == 3 else "Vel"}
n : {params[0]}{"_" if paramsNum == 0 else ""}
m : {params[1]}{"_" if paramsNum == 1 else ""}
Size : {params[2]}{"_" if paramsNum == 2 else ""}
Enter to Continue"""
    
    menu = menu.split("\n")
    for string, x in zip(menu, range(10 + 15 * (len(menu)-1), 9, -15)):
        glRasterPos2i(10, x)
        for c in string:
            glutBitmapCharacter(GLUT_BITMAP_8_BY_13, ord(c))

def activeOverlay():
    menu = f"""StepsPerFrame: {stepsPerFrame}; Inc with 1; Dec with 2
Smooth Colour: {smooth}; Toggle with 3
Velocities: {velocities}; Toggle with 4
Accelerations: {accelerations}; Toggle with 5"""
    
    menu = menu.split("\n")
    for string, x in zip(menu, range(10 + 15 * (len(menu)-1), 9, -15)):
        glRasterPos2i(10, x)
        for c in string:
            glutBitmapCharacter(GLUT_BITMAP_8_BY_13, ord(c))

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
    
    if setup:
        setupOverlay()
    else:
        activeOverlay()
            
    glEnable(GL_DEPTH_TEST)

    glMatrixMode(GL_MODELVIEW)
    glPopMatrix()
    glMatrixMode(GL_PROJECTION)
    glPopMatrix()

def draw():
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glLoadIdentity()
    refresh2d(width, height)

    if setupStage > 0:
        v.interact()
        v.draw()

    overlay()

    glutSwapBuffers()

if __name__ == "__main__":    
    glutInit()
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH)
    glutInitWindowSize(width, height)
    glutInitWindowPosition(100, 100)
    glutCreateWindow(b"Title")
    glutDisplayFunc(draw)
    glutIdleFunc(draw)
    glutKeyboardFunc(handleKeyInput)
    glutMainLoop()