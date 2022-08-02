"""
	Diff.py
	
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
				
if __name__ == "__main__":
		
	def g(x):
		return x*sin(x)

	def dgdx(x):
		return sin(x) + x*cos(x)
	
	def d2gdx2(x):
		return 2.0*cos(x) - x*sin(x)
			
	print("Numerical First Order Derivative:")		
	dg = Diff_1st (g, 0.01)
	print("Forward:", dg.Forward(0.25))
	print("Backward:", dg.Backward(0.25))
	print("Central:", dg.Central(0.25))
	print("Exact: ", dgdx(0.25),"\n")
	
	print("Numerical Second Order Derivative:")	
	d2g = Diff_2nd (g, 0.01)
	print("Forward:", d2g.Forward(0.25))
	print("Backward:", d2g.Backward(0.25))
	print("Central:", d2g.Central(0.25))
	print("Exact: ", d2gdx2(0.25),"\n")
	
