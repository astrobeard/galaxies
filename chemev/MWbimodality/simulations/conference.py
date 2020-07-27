r""" 
Runs a simulation of the inside-out disk model to present at the SDSS 
Gotham 2020 conference. 

ARGV 
----
1) 		The name of the output 
2) 		The number of star particles per zone per timestep 
""" 

import tracers 
import gas_disks 
import common 
import numpy as np 
import math as m 
import vice 
from vice.yields.presets import JW20 
import sys 
import os 

TIME_BINS = np.linspace(0, 12.8, 41).tolist() 
RAD_BINS = np.linspace(0, 30, 121).tolist() 
TIME_SWITCH = 1.0 # time of linear to exponential switch in Gyr 
TSTAR_1 = 0.08 # tau_star normalization in linear phase 
TSTAR_2 = 2. # tau_star normalization in exponential phase 
ZONE_WIDTH = 0.25 # width of each zone in kpc 
RSCALE = 3 # scale radius of this disk model 
DT = 0.01 # The timestep size in Gyr 
TSTAR_NORM = 0.05 
SWITCH = 1.0 


def tau_in(rgal): 
	r""" 
	Star formation timescale in Gyr as a function of galactocentric radius in 
	kpc. 
	""" 
	return 3 + (rgal + 1e-12) / 2.5 


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


class sfe: 

	r""" 
	SFE timescale as a function of time at a given galactocentric radius. 
	""" 

	def __init__(self, rgal): 
		self._norm = tau_star(rgal,  norm = TSTAR_NORM) 
		self._switch = TIME_SWITCH 
		# self._tstar1 = tau_star(rgal, norm = TSTAR_1) 
		# self._tstar2 = tau_star(rgal, norm = TSTAR_2) 
		# self._rgal = rgal 

	def __call__(self, time): 
		return self._norm 
		# if time < self._switch: 
		# 	return self._norm  
		# else: 
		# 	return self._norm + (time - self._switch) / (
		# 		12.8 - self._switch) 


class star_formation_history: 

	r""" 
	The functional form of the star formation history under this model 
	""" 
	def __init__(self, rgal): 
		self._rgal = rgal 
		# self._norm = (2 * m.pi * rgal * m.exp(-rgal / RSCALE) * ZONE_WIDTH / 
		# 	TIME_SWITCH * (
		# 		# TIME_SWITCH / 2. + TSTAR_1 / TSTAR_2 * tau_in(rgal) * (
		# 		TIME_SWITCH / 2. +  tau_in(rgal) * (
		# 			m.exp(-TIME_SWITCH / tau_in(rgal)) - 
		# 			m.exp(-TIME_BINS[-1] / tau_in(rgal)) 
		# 		) 
		# 	)**(-1) 
		# ) 
		self._norm = (2 * m.pi * rgal * m.exp(-rgal / RSCALE) * ZONE_WIDTH / 
			TIME_SWITCH * (
				TIME_SWITCH / 2 + 
				tau_in(rgal) * (
					1 - m.exp(-(12.8 - TIME_SWITCH) / tau_in(rgal))
				) 
			)**(-1) 
		) 

	def __call__(self, time): 
		if time < TIME_SWITCH: 
			# return self._norm 
			return self._norm * time 
		else: 
			# return self._norm * TIME_SWITCH * TSTAR_1 / TSTAR_2 * m.exp(
			return self._norm * TIME_SWITCH * m.exp(
				-(time - TIME_SWITCH) / tau_in(self._rgal)) 


class infall_history: 

	r""" 
	The functional form of the infall history in Msun/yr at a given 
	galactocentric radius 
	""" 
	def __init__(self, rgal): 
		tau_dep = sfe(rgal)(12.8) / (1 + eta(rgal) - 0.4) 
		tau_in_ = tau_in(rgal) 
		self._norm = 2 * m.pi * rgal * m.exp(-rgal / RSCALE) * ZONE_WIDTH * (
			1 - m.exp(-SWITCH / tau_dep) + 
			(tau_in_ / abs(tau_dep - tau_in_)) * (1 - m.exp(-(12.8 - SWITCH) / 
				tau_dep))
		)**(-1) / tau_dep 
		self._tau_in = tau_in_ 
		# self._norm = 2 * m.pi * rgal * m.exp(-rgal / RSCALE) * ZONE_WIDTH * (
		# 	1 - m.exp(-12.8 / tau_dep) 
		# ) / tau_dep 

	def __call__(self, time): 
		# return self._norm 
		if time < SWITCH: 
			return self._norm 
		else: 
			return self._norm * m.exp(-(time - SWITCH) / self._tau_in) 


class diskmodel(vice.multizone): 

	def __init__(self): 
		super().__init__(
			name = sys.argv[1], 
			n_zones = len(RAD_BINS) - 1, 
			n_stars = int(sys.argv[2]), 
			verbose = True, 
			simple = False) 
		self.migration.stars = tracers.UWhydro(TIME_BINS, RAD_BINS, 
			n_stars = self.n_stars, 
			filename = "%s_extra_tracer_data.out" % (self.name)) 
		for i in range(self.n_zones): 
			# self.zones[i].func = infall_history(ZONE_WIDTH * (i + 0.5)) 
			# self.zones[i].mode = "ifr" 
			self.zones[i].func = star_formation_history(ZONE_WIDTH * (i + 0.5)) 
			self.zones[i].mode = "sfr" 
			self.zones[i].bins = np.linspace(-3, 1, 401) 
			self.zones[i].elements = ["fe", "o"] 
			self.zones[i].dt = DT 
			self.zones[i].schmidt = True 
			self.zones[i].Mg0 = 0 
			if i > 61: 
				self.zones[i].func = lambda t: 0 
				# self.zones[i].tau_star = float("inf") 
				self.zones[i].tau_star = 1e4 
				self.zones[i].eta = 100 
				for j in self.zones[i].elements: 
					self.zones[i].entrainment.agb[j] = 0 
					self.zones[i].entrainment.ccsne[j] = 0 
					self.zones[i].entrainment.sneia[j] = 0 
			else: 
				self.zones[i].tau_star = sfe(ZONE_WIDTH * (i + 0.5)) 
				self.zones[i].eta = eta(ZONE_WIDTH * (i + 0.5)) 
				# self.zones[i].Mg0 = 0 
				# self.zones[i].Mg0 = (
				# 	1e9 * self.zones[i].func(0) * 
				# 	self.zones[i].tau_star(0) / 
				# 	(1 + self.zones[i].eta - 0.4) * 
				# 	self.zones[i].MgSchmidt**self.zones[i].schmidt_index 
				# )**(1 / (1 + self.zones[i].schmidt_index)) 

	def run(self): 
		super().run(np.linspace(0, 12.8, 257), overwrite = True) 
		self.migration.stars.close_file() 
		# pass 


def tau_star(rgal, norm = 2): 
	r""" 
	The mathematical relation between galactocentric radius in kpc and the 
	SFE timescale :math:`\tau_\star`. 

	Parameters 
	----------
	rgal : real number 
		Galactocentric radius in kpc 
	norm : real number [default : 2] 
		The value of :math:`\tau_\star` at r = 0. 

	Returns 
	-------
	t : real number 
		:math:`\tau_\star` at the specified radius 

	Notes 
	----- 
	The mathematical relation describing :math:`\tau_\star`: 

	.. math:: \tau_\star = A e^{-r/2r_\text{s}} 

	where :math:`r_\text{s}` is the scale radius, :math:`A` is the norm, 
	and :math:`r` is the galactocentric radius. 
	""" 
	return norm * m.exp(rgal / (2 * RSCALE)) 


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


if __name__ == "__main__": 
	for i in RAD_BINS[:60]: 
		print("R = %.2f kpc ; tau_in = %.2f" % (
			i + ZONE_WIDTH / 2, tau_in(i + ZONE_WIDTH / 2))) 
	diskmodel().run() 


