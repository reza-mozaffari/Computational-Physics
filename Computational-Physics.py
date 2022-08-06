"""
	Computational-Physics.py
	
	M. Reza Mozaffari
	Physics Group, University of Qom
"""

from math import *

class Diff_1st ():
    def __init__ (self, f, h=0.01):
        self.f = f
        self.h = h

    # 2-point O(h)
    def Forward (self, x):
        return (self.f(x + self.h) - self.f(x))/self.h

    # 2-point O(h)
    def Backward (self, x):
        return (self.f(x) - self.f(x - self.h))/self.h

    # 3-point O(h^2)
    def Central (self, x):
        return 0.5*(self.f(x + self.h) - self.f(x - self.h))/self.h
        
class Diff_2nd ():
    def __init__ (self, f, h=0.01):
        self.f = f
        self.h = h

    # 3-point O(h)
    def Forward (self, x):
        return (self.f(x + 2.0*self.h) - 2.0*self.f(x + self.h) + self.f(x))/self.h**2

    # 3-point O(h)
    def Backward (self, x):
        return (self.f(x - 2.0*self.h) - 2.0*self.f(x - self.h) + self.f(x))/self.h**2

    # 3-point O(h^2)
    def Central (self, x):
        return (self.f(x + self.h) - 2.0*self.f(x) + self.f(x - self.h))/self.h**2
        
class Trapezoidal ():
    def __init__ (self, f, a, b, n):
        self.f = f
        self.a = a
        self.b = b
        self.n = n
        self.h = (self.b - self.a)/self.n

    def Integral (self):
        s = 0.0
        for i in range(0, self.n):
            s += (self.f(self.a + i*self.h) + self.f(self.a + (i+1)*self.h))*self.h/2.0
        return s
        
class Simpson ():
    def __init__ (self, f, a, b, n):
        self.f = f
        self.a = a
        self.b = b
        if n%2 == 0:
            self.n = n
        else:
            self.n = n + 1
        self.h = (self.b - self.a)/self.n

    def Integral (self):
        s = 0.0
        for i in range(0, self.n, 2):
            s += ( self.f(self.a + i*self.h) +  4.0*self.f(self.a + (i+1)*self.h) +  self.f(self.a + (i+2)*self.h))*self.h/3.0
        return s
        
class GaussianQuadrature ():
    def __init__ (self, f, a, b, n, w, s):
        self.f = f
        self.a = a
        self.b = b
        self.n = n
        self.h = (self.b - self.a)/self.n
        
        assert len(w) == len(s), "Error in input arguments."
        self.w = []; self.s = []
        for i in range(len(w)):
            self.w.append(w[i])
            self.s.append(s[i])  
        
    def Integral (self):
        s = 0.0
        for i in range(self.n):
            for j in range(len(self.w)):
                s += self.w[j]*self.f(self.a + 0.5*(2*i+1)*self.h + 0.5*self.h*self.s[j])*self.h/2.0
        return s
        
class Bisection():
    def __init__ (self, f, a, b, tol=1.0e-3):
        self.f = f
        self.a = a
        self.b = b
        self.tol = tol
        
        assert self.f(a)*self.f(b)<0, "Error in input arguments."
        
    def FindingRoot(self):
        a, b = self.a, self.b
        c = 0.5*(a + b)
        while abs(self.f(c)) > self.tol:
            if self.f(a)*self.f(c)>0:
                a = c
            elif self.f(a)*self.f(c)<0: 
                b = c
            c = 0.5*(a + b)
        
        return c
        
class NewtonRaphson():
    def __init__ (self, f, fp, x0, tol=1.0e-3):
        self.f = f
        self.fp = fp
        self.x0 = x0
        self.tol = tol
                
    def FindingRoot(self):
        x0 = self.x0
        while abs(self.f(x0)) > self.tol:
            x0 = x0 - self.f(x0)/self.fp(x0)
        return x0


