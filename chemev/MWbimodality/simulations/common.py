r""" 
Routines common to the multizone simulations in written here. 
""" 

__all__ = ["TIME_BINS", "RAD_BINS", "ZONE_WIDTH", "eta", "tau_star", 
	"harmonic_timescale", "depletion_time"]   
import vice 
from vice.yields.presets import my_yields 
import numpy as np 
import math as m 

TIME_BINS = np.linspace(0, 12.8, 41).tolist() 
RAD_BINS = np.linspace(0, 30, 121).tolist() 
ZONE_WIDTH = 0.25 


def eta(rgal, corrective = 0): 
	r""" 
	The mathematical relation between galactocentric radius in kpc and the 
	mass loading factor. 

	Parameters 
	----------
	rgal : real number 
		Galactocentric radius in kpc 
	corrective : real number [default : 0] 
		The corrective term to account for tau_star / tau_sfh 

	Returns 
	-------
	eta : real number 
		The mass loading factor at that radius. 
	""" 
	return vice.yields.ccsne.settings['o'] / vice.solar_z['o'] * (
		10**(0.06 * (rgal - 4) - 0.3)) - 0.6 + corrective


def tau_star(rgal, norm = 2, scale = 3): 
	r""" 
	The mathematical relation between galactocentric radius in kpc and the 
	SFE timescale :math:`\tau_\star`. 

	Parameters 
	----------
	rgal : real number 
		Galactocentric radius in kpc 
	norm : real number [default : 2] 
		The value of :math:`\tau_\star` at r = 0. 
	scale : real number [default : 3] 
		The scale radius of the gas disk 

	Returns 
	-------
	t : real number 
		:math:`\tau_\star` at the specified radius 

	Notes 
	----- 
	The mathematical relation describing :math:`\tau_\star`: 

	.. math:: \tau_\star = A e^{-r/2r_\text{s}} 

	where :math:`r_\text{s}` is the scale radius, :math:`A` is the norm, and 
	:math:`r` is the galactocentric radius. 
	""" 
	# return norm * m.exp(rgal / (2 * scale)) 
	return norm * m.exp(rgal / (5 * scale)) 


def harmonic_timescale(t1, t2): 
	r""" 
	Calculate a harmonic timescale. 

	Parameters 
	----------
	t1 : real number 
		The first timescale. 
	t2 : real number 
		The second timescale. 

	Returns 
	-------
	tau : real number 
		The harmonic timescale defined by: 

		.. math:: \tau \equiv \left(t1^{-1} - t2^{-1}\right)^{-1} 
	""" 
	if t1 and t2: 
		return (1 / t1 - 1 / t2)**(-1) 
	else: 
		return 0 


def depletion_time(tau_star, eta = 2.5, r = 0.4): 
	r""" 
	Calculate a depletion time 

	Parameters 
	----------
	tau_star : real number 
		The star formation efficiency timescale 
	eta : real number [default : 2.5] 
		The mass loading factor 
	r : real number [default : 0.4] 
		The recycling parameter 

	Returns 
	-------
	tau_dep : real number 
		The depletion timescale defined by: 

		.. math:: \tau_\text{dep} \equiv \frac{\tau_\star}{1 + \eta - r} 
	""" 
	return tau_star / (1 + eta - r) 
