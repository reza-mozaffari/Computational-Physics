"""
	Trapezoidal.py
	
	M. Reza Mozaffari
	Physics Group, University of Qom
"""
from math import *

class Trapezoidal ():
	def __init__ (self, f, a, b, n):
		self.f = f
		self.a = a
		self.b = b
		self.n = n
		self.h = (b - a)/n
	
	def Integral (self):
		s = 0.0
		for i in range(0, self.n):
			s += ( self.f(self.a + i*self.h) + self.f(self.a + (i+1)*self.h) )*self.h/2.0
		return s

if __name__ == "__main__":
		
	def g(x):
		return x*sin(x)

	trap = Trapezoidal (g, 0.0, pi, 250)
	print(trap.Integral())
