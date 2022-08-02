"""
	Simpson.py
	
	M. Reza Mozaffari
	Physics Group, University of Qom
"""
from math import *

class Simpson ():
	def __init__ (self, f, a, b, n):
		self.f = f
		self.a = a
		self.b = b
		if n%2 == 0:
			self.n = n
		else:
			self.n = n + 1
		self.h = (b - a)/n
	
	def Integral (self):
		s = 0.0
		for i in range(0, self.n, 2):
			s += ( self.f(self.a + i*self.h) +  4.0*self.f(self.a + (i+1)*self.h) +  self.f(self.a + (i+2)*self.h))*self.h/3.0
		return s

if __name__ == "__main__":
		
	def g(x):
		return x*sin(x)

	simp = Simpson (g, 0.0, pi, 100)
	print(simp.Integral(), pi)
