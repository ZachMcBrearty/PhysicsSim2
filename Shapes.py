from math import sin, cos, radians, atan

from OpenGL.GL import glBegin, glEnd, glColor3fv, glVertex2fv, glLineWidth, \
    GL_POLYGON, GL_QUADS, GL_TRIANGLES, GL_LINES

class Vector:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    
    def rotateby(self, angle):
       angle = radians(angle)
       x = self.x*cos(angle) - self.y*sin(angle)
       y = self.x*sin(angle) + self.y*cos(angle)
        
       return Vector(x, y)
    
    @classmethod
    def frompolar(cls, mag, angle):
        angle = radians(angle)
        return cls(cos(angle), sin(angle)) * mag
    
    def __str__(self):
        return f"({self.x}, {self.y})"
    __repr__=__str__
    def __add__(self, other):
        if isinstance(other, Vector):
            x = self.x + other.x
            y = self.y + other.y
            return Vector(x, y)
        else:
            return NotImplemented
    
    def __sub__(self, other):
        if isinstance(other, Vector):
            x = self.x - other.x
            y = self.y - other.y
            return Vector(x, y)
        else:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            x = self.x * other
            y = self.y * other
            return Vector(x, y)
        else:
            return NotImplemented
    def __rmul__(self, other):
        if isinstance(other, (int, float)):
            x = self.x * other
            y = self.y * other
            return Vector(x, y)
        else:
            return NotImplemented
        
    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            if other == 0:
                return Vector(0, 0)
            else:
                x = self.x / other
                y = self.y / other
                return Vector(x, y)
        else:
            return NotImplemented
        
    def __rtruediv__(self, other):
        if isinstance(other, (int, float)):
            if self.x == 0:
                return Vector(0, 0)
            else:
                x = self.x / other
                y = self.y / other
                return Vector(x, y)
        else:
            return NotImplemented
        
    def __abs__(self):
        return (self.x**2 + self.y**2)**0.5
    
    def normalise(self):
        return self / abs(self)

class Polygon:
    def __init__(self, x, y, vertices, fillcolour=None, strokecolour=None, strokewidth=4):
        '''baseclass all shapes,
        if the colours are None then they are not drawn'''
        self.pos = Vector(x, y)
        self.vertices = vertices
        self.fillcolour = fillcolour
        self.strokecolour = strokecolour
        self.shapetype = GL_POLYGON
        self.strokewidth = strokewidth
        
    def draw(self):
        if self.fillcolour is not None:
            glBegin(self.shapetype)
            glColor3fv(self.fillcolour)
            for v in self.vertices:
                glVertex2fv(v)
            glEnd()
        
        if self.strokecolour is not None:
            glLineWidth(self.strokewidth)
            glBegin(GL_LINES)
            glColor3fv(self.strokecolour)
            # each glVertex2fv _pair_ creates a single line
            for i in range(len(self.vertices)-1):
                glVertex2fv(self.vertices[i])
                glVertex2fv(self.vertices[i+1])
            # pair the first and last to finish the shape
            glVertex2fv(self.vertices[-1])
            glVertex2fv(self.vertices[0])
            glEnd()
            
    def move(self, x, y):
        self.pos += Vector(x, y)
        for v in self.vertices:
            v[0] += x
            v[1] += y
            
    def dist(self, other):    
        return self.pos - other.pos
        
class Circle(Polygon):
    def __init__(self, x, y, rad, res, 
                 fillcolour=None, strokecolour=None, strokewidth=4):
        self.pos = Vector(x, y)
        self.rad = rad
        self.res = res
        
        self.vertices = [[self.pos.x + self.rad * sin(radians(angle)), 
                          self.pos.y + self.rad * cos(radians(angle))]
                         for angle in range(-180, 181, self.res)]
        
        self.fillcolour = fillcolour
        self.strokecolour = strokecolour
        self.shapetype = GL_POLYGON
        self.strokewidth = strokewidth
        
if __name__ == "__main__":
    a = Vector(1, 2)
    b = Vector(3, 4)
    print("a", a)
    print("b", b)
    print("a+b", a+b)
    
    a += b
    
    print("a += b", a)