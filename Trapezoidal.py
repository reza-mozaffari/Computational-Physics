"""
	Trapezoidal.py
	
	M. Reza Mozaffari
	Physics Group, University of Qom
"""

class Trapezoidal ():
	def __init__ (self, f, a, b, n):
		self.f = f
		self.a = a
		self.b = b
		self.n = n
		self.h = (b - a)/n
	
	def Integral (self):
		s = 0.5*(self.f(self.a) + self.f(self.b))*self.h
		for i in range(1,self.n):
			s = s + self.f(self.a + i*self.h)*self.h

		return s

if __name__ == "__main__":
		
	def g(x):
		return x

	trap = Trapezoidal (g, 1.0, 2.0, 100)
	print(trap.Integral())
