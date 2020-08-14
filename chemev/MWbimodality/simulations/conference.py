r""" 
Runs a simulation of the inside-out disk model to present at the SDSS 
Gotham 2020 conference. 

ARGV 
----
1) 		The name of the output 
2) 		The number of star particles per zone per timestep 
""" 

# import tracers 
import gas_disks 
import common 
import numpy as np 
import math as m 
import vice 
from vice.yields.presets import JW20  
from vice.toolkit import hydrodisk 
import sys 
import os 

TIME_BINS = np.linspace(0, 12.8, 41).tolist() 
RAD_BINS = np.linspace(0, 30, 121).tolist() 
TIME_SWITCH = 2. # time of linear to exponential switch in Gyr 
TSTAR_1 = 0.04 # tau_star normalization in linear phase 
TSTAR_2 = 2. # tau_star normalization in exponential phase 
ZONE_WIDTH = 0.25 # width of each zone in kpc 
RSCALE = 3 # scale radius of this disk model 
DT = 0.01 # The timestep size in Gyr 
TSTAR_NORM = 0.2 
ALPHA = 0.1 


FMOL = 0.35 # the ratio of molecular to total gas content in the disk 
RS_MOL = 3. # scale radius of molecular disk in kpc 
RS_HI = 6. # scale radius of the HI disk in kpc 
R_SF = 15.5 # The radius of the star forming disk in kpc 
R_CMZ = 0.5 # radius of the central molecular zone in kpc 
R_HI = 20 # radius of the HI disk in kpc 
TAU_STAR_MOL = 2.45 # depletion timescale of molecular gas 
TSWITCH = 2. # linear to exponential switch in Gyr 
REFF = 6. # The present-day effective radius of the disk galaxy 
# MGSCHMIDT0 = 1.e6 

SANCHEZ_TAU_SFH_DATA = np.genfromtxt("sanchez_tau_sfh.dat").tolist() 


def interpolate(x1, x2, y1, y2, x): 
	r""" 
	A simple 1-dimensional interpolation routine derived from the point-slope 
	form of a line. Interpolate the y-value of a point x given two points on 
	the xy-axis (x1, y1) and (x2, y2). 
	""" 
	return (y2 - y1) / (x2 - x1) * (x - x1) + y1 


def to_Re(rgal, Re = 6): 
	r""" 
	Convert a galactocentric radius in kpc to a galactocentric radius in units 
	of the effective radius Re. 

	rgal : galactocentric radius in kpc 
	Re : effective radius in kpc [default : 6] 
	""" 
	return rgal / Re 


def get_bin_number(val, bins): 
	for i in range(len(bins) - 1): 
		if bins[i] <= val <= bins[i + 1]: return i 
	return -1 


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
		if rgal <= R_CMZ: 
			self._tstar = TAU_STAR_MOL 
		elif rgal <= R_SF: 
			self._tstar = TAU_STAR_MOL * (1 + sfe.amp_ratio() * 
				m.exp(-(rgal - R_CMZ) / RS_HI) / 
				m.exp(-rgal / RS_MOL) 
			) 
		else: 
			self._tstar = 1.e6 

	def __call__(self, time): 
		return self._tstar 

	@staticmethod 
	def amp_ratio(): 
		r""" 
		The amplitude ratio of the HI to molecular gas disks 
		""" 
		return (1 / FMOL - 1) * (
			RS_MOL**2 * (1 - m.exp(-R_SF / RS_MOL)) - RS_MOL * R_SF * m.exp(
				-R_SF / RS_MOL) 
		) / (
			R_CMZ * RS_HI - R_HI * RS_HI * m.exp(-(R_HI - R_CMZ) / RS_HI) + 
			RS_HI**2 * (1 - m.exp(-(R_HI - R_CMZ) / RS_HI)) 
		) 


# class star_formation_history: 

# 	r""" 
# 	The functional form of the star formation history under this model 
# 	""" 
# 	def __init__(self, rgal): 
# 		self._rgal = rgal 
# 		self._norm = (2 * m.pi * rgal * m.exp(-rgal / RSCALE) * ZONE_WIDTH / 
# 			TIME_SWITCH * (
# 				# TIME_SWITCH / 2. + TSTAR_1 / TSTAR_2 * tau_in(rgal) * (
# 				TIME_SWITCH / 2. +  tau_in(rgal) * (
# 					m.exp(-TIME_SWITCH / tau_in(rgal)) - 
# 					m.exp(-TIME_BINS[-1] / tau_in(rgal)) 
# 				) 
# 			)**(-1) 
# 		) 

# 	def __call__(self, time): 
# 		if time < TIME_SWITCH: 
# 			# return self._norm 
# 			return self._norm * TIME_SWITCH  
# 		else: 
# 			# return self._norm * TIME_SWITCH * TSTAR_1 / TSTAR_2 * m.exp(
# 			return self._norm * TIME_SWITCH * m.exp(
# 				-(time - TIME_SWITCH) / tau_in(self._rgal)) 


class star_formation_history: 

	def __init__(self, rgal): 
		self._timescale = star_formation_history.tau_sfh(rgal) 
		self._norm = self.norm(rgal) 

	def __call__(self, time): 
		# linear-exponential form 
		# return self._norm * time * m.exp(-time / self._timescale) 

		# linear-then-exponential form 
		if time <= TSWITCH: 
			return self._norm * time 
		else: 
			return self._norm * TSWITCH * m.exp(-(time - TSWITCH) / 
				self._timescale) 

	@property 
	def timescale(self): 
		return self._timescale 

	@staticmethod 
	def tau_sfh(rgal): 
		r""" 
		Star formation history timescale as a function of galactocentric 
		radius. 

		Parameters 
		----------
		rgal : real number  
			Galactocentric radius in kpc 

		Returns 
		-------
		tau_sfh : real number  
			The e-folding timescale of the star formation history at that 
			radius in Gyr. 
		""" 
		# The simple linear-with-radius form used previously 
		# return 3 + (rgal + 1e-12) / 2.5 

		# The Sanchez (2020) timescales 
		Re = [i[0] for i in SANCHEZ_TAU_SFH_DATA] 
		tau = [i[1] for i in SANCHEZ_TAU_SFH_DATA] 
		rgal = to_Re(rgal) 
		bin_ = get_bin_number(rgal, Re) 
		if bin_ != -1: 
			return interpolate(Re[bin_], Re[bin_ + 1], tau[bin_], tau[bin_ + 1], 
				rgal) 
		else: 
			return interpolate(Re[-2], Re[-1], tau[-2], tau[-1], rgal) 


	def norm(self, rgal): 
		r""" 
		The normalization of the star formation history as a function of 
		galactocentric radius. 

		Parameters 
		----------
		rgal : real number 
			Galactocentric radius in kpc 

		Returns 
		-------
		norm : real number 
			The normalization of the star formation history in 
			:math:`M_\odot/yr`. 
		""" 

		# norm for linear-exponential evolution 
		# return self.timescale**(-2) * (
		# 	1 - (1 + TIME_BINS[-1] / self.timescale) * m.exp(
		# 		-TIME_BINS[-1] / self.timescale
		# 	) 
		# )**(-1) * rgal * m.exp(-rgal / RS_MOL) 

		# norm for linear-then-exponential evolution 
		return rgal * m.exp(-rgal / RS_MOL) * ZONE_WIDTH * (
			0.5 * TSWITCH**2 + 
			TSWITCH * star_formation_history.tau_sfh(rgal) * (
				1 - m.exp(-(TIME_BINS[-1] - TSWITCH) / 
					star_formation_history.tau_sfh(rgal)
				) 
			) 
		)**(-1) 


class fiducial_sfh: 

	def __init__(self, rgal): 
		self._tau_rise = 2 
		self._timescale = star_formation_history.tau_sfh(rgal) 
		self._norm = rgal * m.exp(-rgal / RSCALE) * (
			self._timescale * (1 - m.exp(-12.8 / self._timescale)) - 
			self._tau_rise * self._timescale / (
				self._tau_rise + self._timescale) * 
			(1 - m.exp(-12.8 * (self._tau_rise + self._timescale) / (
				self._tau_rise * self._timescale))) 
		)**(-1) * 3 

	def __call__(self, time): 
		return self._norm * m.exp(-time / self._timescale) * (1 - 
			m.exp(-time / self._tau_rise)) 

	@property 
	def timescale(self): 
		return self._timescale 


class fiducial_sfh_with_lateburst(fiducial_sfh): 

	def __init__(self, rgal): 
		super().__init__(rgal) 
		self._tmax = 10.8
		self._width = 1 
		self._a = super().__call__(self._tmax) 

	def __call__(self, time): 
		return 0.91 * (super().__call__(time) + self._a * m.exp(
			-(time - self._tmax)**2 / (2 * self._width)**2 
		))


class constant_sfh: 

	r""" 
	Constant star formation history embedded within the molecular plus 
	neutral hydrogen star-forming disk. 
	""" 

	def __init__(self, rgal): 
		self._mdotstar0 = 1. 
		self._rgal = rgal 

	def __call__(self, time): 
		return self._mdotstar0 * 2 * m.pi * self._rgal * m.exp(-self._rgal / 
			RS_MOL) * ZONE_WIDTH 


class constant_gas: 

	def __init__(self, rgal): 
		self._mass = 1.e8 * rgal * m.exp(-rgal / RSCALE) 

	def __call__(self, time): 
		return self._mass 


class infall_history: 

	r""" 
	The functional form of the infall history in Msun/yr at a given 
	galactocentric radius 
	""" 
	def __init__(self, rgal): 
		# self._rgal = rgal 
		# self._tstar = tau_star(rgal, norm = TSTAR_NORM)
		# self._eta = eta(rgal, 
		# 	corrective = self._tstar / tau_in(rgal)) 
		# self._tau_dep = self._tstar / (1 + self._eta - 0.4) 
		# self._norm = 2 * m.pi * rgal * m.exp(-rgal / RSCALE) * ZONE_WIDTH * (
		# 	harmonic_timescale(self._tau_dep, tau_in(rgal)) * (
		# 		m.exp(-12.8 / tau_in(rgal)) - 
		# 		m.exp(-12.8 / self._tau_dep)
		# 	) + ALPHA * self._tau_dep * (
		# 		1 - m.exp(-12.8 / self._tau_dep) 
		# 	)
		# )**(-1) 
		# self._tau_in = tau_in(rgal) 
		self._tau_in = star_formation_history.tau_sfh(rgal) 
		self._tstar = TAU_STAR_MOL 
		self._eta = eta(rgal, corrective = self._tstar / self._tau_in) 
		self._tau_dep = self._tstar / (1 + self._eta - 0.4) 
		self._norm = 1.e-3 * rgal * m.exp(-rgal / RSCALE) * harmonic_timescale(
			self._tau_in, self._tau_dep)**(-1) * (
			m.exp(-12.8 / self._tau_dep) - 
			m.exp(-12.8 / self._tau_in) 
		)**(-1) 


	def __call__(self, time): 
		return self._norm * m.exp(-time / self._tau_in) 

	@property 
	def timescale(self): 
		return self._tau_in 



class diskmigration(hydrodisk.hydrodiskstars): 

	r""" 
	Subclassed hydrodiskstars object to write extra analog star particle data 
	to an output file. 
	""" 

	def __init__(self, radbins, mode = "linear", filename = "stars.out"): 
		super().__init__(radbins, mode = mode) 
		if isinstance(filename, str): 
			self._file = open(filename, 'w') 
			self._file.write("# zone_origin\ttime_origin\tzfinal\n") 
		else: 
			raise TypeError("Filename must be a string. Got: %s" % (
				type(filename))) 

		# Multizone object automatically swaps this to True in setting up 
		# its stellar population zone histories 
		self.write = False 

	def __call__(self, zone, tform, time): 
		if tform == time: 
			super().__call__(zone, tform, time) # reset analog star particle 
			if self.write: 
				if self.analog_index == -1: 
					finalz = 100 
				else: 
					finalz = self.analog_data["zfinal"][self.analog_index] 
				self._file.write("%d\t%.2f\t%.2f\n" % (zone, tform, finalz)) 
			else: pass 
			return zone 
		else: 
			return super().__call__(zone, tform, time) 

	def close_file(self): 
		r""" 
		Closes the output file - should be called after the multizone model 
		simulation runs. 
		""" 
		self._file.close() 

	@property 
	def write(self): 
		r""" 
		Type : bool 

		Whether or not to write out to the extra star particle data output 
		file. For internal use by the vice.multizone object only. 
		""" 
		return self._write 

	@write.setter 
	def write(self, value): 
		if isinstance(value, bool): 
			self._write = value 
		else: 
			raise TypeError("Must be a boolean. Got: %s" % (type(value))) 


class diskmodel(vice.multizone): 

	def __init__(self): 
		super().__init__(
			name = sys.argv[1], 
			n_zones = len(RAD_BINS) - 1, 
			n_stars = int(sys.argv[2]), 
			verbose = True, 
			simple = False) 
		# self.migration.stars = tracers.UWhydro(TIME_BINS, RAD_BINS, 
		# 	n_stars = self.n_stars, 
		# 	filename = "%s_extra_tracer_data.out" % (self.name)) 
		# self.migration.stars = hydrodisk.hydrodiskstars(RAD_BINS) 
		self.migration.stars = diskmigration(RAD_BINS, 
			filename = "%s_extra_tracer_data.out" % (sys.argv[1])) 
		for i in range(self.n_zones): 
			# self.zones[i].func = infall_history(ZONE_WIDTH * (i + 0.5)) 
			# self.zones[i].mode = "ifr" 
			# self.zones[i].func = star_formation_history(ZONE_WIDTH * (i + 0.5)) 
			# self.zones[i].func = constant_sfh(ZONE_WIDTH * (i + 0.5)) 
			# self.zones[i].func = fiducial_sfh(ZONE_WIDTH * (i + 0.5)) 
			self.zones[i].func = fiducial_sfh_with_lateburst(
				ZONE_WIDTH * (i + 0.5)) 
			self.zones[i].mode = "sfr" 
			# self.zones[i].func = constant_gas(ZONE_WIDTH * (i + 0.5)) 
			# self.zones[i].mode = "gas" 
			self.bins = np.linspace(-3, 1, 401) 
			self.zones[i].elements = ["fe", "o"] 
			self.zones[i].dt = DT 
			self.zones[i].Mg0 = 0 
			self.zones[i].schmidt = True 
			# self.zones[i].schmidt = False 
			# if i: 
			# 	self.zones[i].MgSchmidt = MGSCHMIDT0 * (RAD_BINS[i + 1]**2 - 
			# 		RAD_BINS[i]**2) / RAD_BINS[1]**2 
			# else: 
			# 	self.zones[i].MgSchmidt = MGSCHMIDT0 
			self.zones[i].MgSchmidt = 1.e7 * m.pi * (RAD_BINS[i + 1]**2 - 
				RAD_BINS[i]**2) 
			self.zones[i].schmidt_index = 0.85 
			if i > 61: 
				self.zones[i].func = lambda t: 0 
				self.zones[i].tau_star = 1.e6 
				self.zones[i].eta = 100 
				for j in self.zones[i].elements: 
					self.zones[i].entrainment.agb[j] = 0 
					self.zones[i].entrainment.ccsne[j] = 0 
					self.zones[i].entrainment.sneia[j] = 0 
			else: 
				# self.zones[i].tau_star = sfe(ZONE_WIDTH * (i + 0.5)) 
				# self.zones[i].eta = eta(ZONE_WIDTH * (i + 0.5), 
				# 	corrective = (self.zones[i].tau_star(0) / 
				# 		self.zones[i].func.timescale)
				# 	) 
				self.zones[i].tau_star = TAU_STAR_MOL 
				# self.zones[i].tau_star = tau_star(ZONE_WIDTH * (i + 0.5), 
				# 	norm = TAU_STAR_MOL) 
				self.zones[i].eta = eta(ZONE_WIDTH * (i + 0.5), 
					# corrective = self.zones[0].tau_star / 
					# 	self.zones[i].func.timescale 
					corrective = 0 
					) 

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
		print("R = %.2f kpc ; tau_sfh = %.2f" % (
			i + ZONE_WIDTH / 2, 
			star_formation_history.tau_sfh(i + ZONE_WIDTH / 2))) 
	diskmodel().run() 


