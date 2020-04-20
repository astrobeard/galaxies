
__all__ = ["static_exponential", "exponential_decay"]  
import math as m 
import numbers 

class static_exponential: 

	""" 
	A callable object for a gas reservoir within a static exponential disk 
	""" 

	def __init__(self, zone, norm, rad_bins, scale_length): 
		self._mass = norm * m.exp(-zone * (rad_bins[zone + 1] - rad_bins[zone]) / 
			scale_length) * (rad_bins[zone + 1]**2 - 
			rad_bins[zone]**2) / rad_bins[1]**2 

	def __call__(self, time): 
		return self._mass 


class exponential_decay: 

	r""" 
	A callable object representing exponential decay. 
	""" 

	def __init__(self, norm, timescale): 
		self.norm = norm 
		self.timescale = timescale 

	def __call__(self, t): 
		return self.norm * m.exp(-t / self.timescale) 

	@property 
	def norm(self): 
		r""" 
		The normalization of the exponential (i.e. the value at t = 0) 
		""" 
		return self._norm 

	@norm.setter 
	def norm(self, value): 
		if isinstance(value, numbers.Number): 
			if value > 0: 
				self._norm = float(value)  
			else: 
				raise ValueError("Must be positive. Got: %g" % (value)) 
		else: 
			raise TypeError("Must be a numerical value. Got: %s" % (
				type(value))) 

	@property 
	def timescale(self): 
		r""" 
		The decay timescale of the exponential 
		""" 
		return self._timescale 

	@timescale.setter 
	def timescale(self, value): 
		if isinstance(value, numbers.Number): 
			if value > 0: 
				self._timescale = float(value) 
			else: 
				raise ValueError("Must be positive. Got: %g" % (value)) 
		else: 
			raise TypeError("Must be a numerical value. Got: %s" % (
				type(value))) 
			

class linear_exponential(exponential_decay): 

	r""" 
	A callable object representing a linear-exponential evolution 
	""" 
	def __init__(self, norm, timescale): 
		super().__init__(norm, timescale) 

	def __call__(self, time): 
		return self.norm * time * m.exp(-time / self.timescale) 


class linear_then_exponential(exponential_decay): 

	r""" 
	A callable object representing a linear-then-exponential evolution 
	""" 

	def __init__(self, norm, timescale, switch): 
		super().__init__(norm, timescale) 
		self.switch = switch 

	def __call__(self, time): 
		if time <= self.switch: 
			return self.norm * time 
		else: 
			return self.norm * self.switch * m.exp(-(time - self.switch) / 
				self.timescale) 

	@property 
	def switch(self): 
		r""" 
		The time at which the evolution switches from linear growth to 
		exponential decay. 
		""" 
		return self._switch 

	@switch.setter 
	def switch(self, value): 
		if isinstance(value, numbers.Number): 
			if value > 0: 
				self._switch = float(value) 
			else: 
				raise ValueError("Must be positive. Got: %g" % (value)) 
		else: 
			raise TypeError("Must be a numerical value. Got: %s" % (
				type(value))) 



