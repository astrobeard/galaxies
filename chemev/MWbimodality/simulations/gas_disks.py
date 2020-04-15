
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


class linear_then_exponential: 

	r""" 
	A callable object representing a linear-then-exponential evolution 
	""" 

	def __init__(self, slope, delay, timescale): 
		self.slope = slope 
		self.delay = delay 
		self.timescale = timescale 

	def __call__(self, time): 
		if time < self.delay: 
			return self.slope * time 
		else: 
			return (self.slope * self.delay) * m.exp(
				-(time - self.delay) / self.timescale
			) 

	@property 
	def slope(self): 
		r""" 
		The slope of the increase in Gyr^-1 
		""" 
		return self._slope 

	@slope.setter 
	def slope(self, value): 
		if isinstance(value, numbers.Number): 
			if value > 0: 
				self._slope = value 
			else: 
				raise ValueError("Must be positive. Got: %g" % (value)) 
		else: 
			raise TypeError("Must be a numerical value. Got: %s" % (
				type(value))) 

	@property 
	def delay(self): 
		r""" 
		The amount of time the evolution increases linearly in Gyr 
		""" 
		return self._delay 

	@delay.setter 
	def delay(self, value): 
		if isinstance(value, numbers.Number): 
			if value > 0: 
				self._delay = value 
			else: 
				raise ValueError("Must be positive. Got: %g" % (value)) 
		else: 
			raise TypeError("Must be a numerical value. Got: %s" % (
				type(value))) 

	@property 
	def timescale(self): 
		r""" 
		The e-folding timescale of the exponential decay. 
		""" 
		return self._timescale 

	@timescale.setter 
	def timescale(self, value): 
		if isinstance(value, numbers.Number): 
			if value > 0: 
				self._timescale = value 
			else: 
				raise ValueError("Must be positive. Got: %g" % (value)) 
		else: 
			raise TypeError("Must be a numerical value. Got: %s" % (
				type(value))) 




