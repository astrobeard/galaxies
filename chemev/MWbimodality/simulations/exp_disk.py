""" 
Runs a simulation of an exponential disk with VICE's multizone object 
assuming constant SFR with tracer particles tuned to given output data 

ARGV 
==== 
1)	The relative path to the file containing star particle data 
2) 	The amount of time the hydro simulation ran 
3)	The column number of that file containing star particle formation times 
4) 	The column number of that file containing formation radii 
5) 	The column number of that file containing final radii 
6) 	The scale length of the disk in kpc 
""" 

import numpy as np 
import math as m 
import vice 
import sys 
import os 

TIME_BINS = np.linspace(0, 13.8, 21).tolist() 
RAD_BINS = np.linspace(0, 30, 121).tolist() 
ZONE_WIDTH = 0.25 

def particle_data(filename = "../data/UWhydro_particles.dat"): 
	""" 
	Reads in the particle data from the specified file 
	""" 
	data = np.genfromtxt(filename).tolist() 
	tform = [row[int(sys.argv[3])] for row in data] 
	Rform = [row[int(sys.argv[4])] for row in data] 
	Rfinal = [row[int(sys.argv[5])] for row in data] 
	return [list(i) for i in zip(tform, Rform, Rfinal)] 


def get_bin_number(bins, val): 
	""" 
	Gets the bin number of a given value in a given binspace. 

	bins :: the array containing the bin edges. Assumed to be sorted. 
	val :: the number whose bin number needs to be found 

	Returns -1 if the value does not lie within the binspace 
	""" 
	for i in range(len(bins) - 1): 
		if bins[i] <= val <= bins[i + 1]: return i 
	return -1 


def interpolate(x1, x2, y1, y2, x): 
	try: 
		return (y2 - y1) / (x2 - x1) * (x - x1) + y1 
	except ZeroDivisionError: 
		return y2 

class tracer_settings(object): 

	def __init__(self, time_bins, rad_bins): 
		self._time_bins = time_bins[:] 
		self._rad_bins = rad_bins[:] 
		self._zones = self._analyze_radii(particle_data(sys.argv[1])) 

	def __call__(self, zone, time): 
		tbin = get_bin_number(self._time_bins, time) 
		final = self._zones[zone][tbin][np.random.randint(len(
			self._zones[zone][tbin]))] + np.random.random() 
		init = zone + np.random.random() 
		def zones(t): 
			if t < time: 
				return 0 
			elif t == time: 
				return zone 
			else: 
				return int(interpolate(time, self._time_bins[-1], init, final, 
					t)) 
		return zones 

	def _analyze_radii(self, data): 
		zones = (len(self._rad_bins) - 1) * [None] 
		for i in range(len(zones)): 
			zones[i] = len(self._time_bins) * [None] 
			for j in range(len(zones[i])): 
				zones[i][j] = [] 
		print("Analyzing particle data....") 
		for i in range(len(data)): 
			tbin = get_bin_number(self._time_bins, data[i][0]) 
			rbin = get_bin_number(self._rad_bins, data[i][1]) 
			zones[rbin][tbin].append(get_bin_number(self._rad_bins, data[i][2]))  
			sys.stdout.write("Progress: %.2f%%\r" % (100 * (i + 1) / len(data))) 
			sys.stdout.flush() 
		sys.stdout.write("\n") 
		for i in range(len(zones)): 
			for j in range(len(zones[i])): 
				if len(zones[i][j]) == 0: 
					zones[i][j].append(i) 
				else: 
					continue 
		return zones 

class gas_reservoir(object): 

	""" 
	The callable object for the gas reservoir in a given zone 
	""" 
	def __init__(self, zone): 
		self._mass = 6.0e9 * m.exp(-zone * ZONE_WIDTH / float(sys.argv[6])) * (
			RAD_BINS[zone + 1]**2 - RAD_BINS[zone]**2) / RAD_BINS[1]**2 

	def __call__(self, time): 
		return self._mass 


def run_simulation(): 
	from vice.yields.presets import my_yields 
	mz = vice.multizone(name = "test", n_zones = len(RAD_BINS) - 1, 
		n_tracers = 2, verbose = True, simple = False) 
	mz.migration.stars = tracer_settings(TIME_BINS, RAD_BINS) 
	for i in range(mz.n_zones): 
		mz.zones[i].mode = "gas" 
		mz.zones[i].func = gas_reservoir(i) 
		mz.zones[i].bins = np.linspace(-3, 1, 401) 
		mz.zones[i].elements = ["mg", "fe", "o"] 
		mz.zones[i].eta = 0.2 * (i + 1) 
	print("Running....") 
	mz.run(np.linspace(0, 13.8, 1381), overwrite = True) 

if __name__ == "__main__": 
	run_simulation() 

