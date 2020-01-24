
__all__ = ["static_exponential"] 
import math as m 

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



