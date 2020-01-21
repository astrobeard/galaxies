""" 
Runs a simulation of an exponential disk with VICE's multizone object 
assuming constant SFR with tracer particles tuned to given output data 

ARGV 
==== 
1) 	The scale length of the disk in kpc 
""" 

from tracers import * 
import numpy as np 
import math as m 
import vice 
import sys 
import os 

TIME_BINS = np.linspace(0, 13.8, 21).tolist() 
RAD_BINS = np.linspace(0, 30, 121).tolist() 
ZONE_WIDTH = 0.25 

class gas_reservoir(object): 

	""" 
	The callable object for the gas reservoir in a given zone 
	""" 
	def __init__(self, zone): 
		self._mass = 6.0e9 * m.exp(-zone * ZONE_WIDTH / float(sys.argv[1])) * (
			RAD_BINS[zone + 1]**2 - RAD_BINS[zone]**2) / RAD_BINS[1]**2 

	def __call__(self, time): 
		return self._mass 


def run_simulation(): 
	from vice.yields.presets import my_yields 
	mz = vice.multizone(name = "outward_gas", n_zones = len(RAD_BINS) - 1, 
		n_tracers = 2, verbose = True, simple = False) 
	mz.migration.stars = UWhydro(TIME_BINS, RAD_BINS) 
	for i in range(mz.n_zones): 
		mz.zones[i].mode = "gas" 
		mz.zones[i].func = gas_reservoir(i) 
		mz.zones[i].bins = np.linspace(-3, 1, 401) 
		mz.zones[i].elements = ["mg", "fe", "o"] 
		mz.zones[i].eta = 0.2 * (i + 1) 
		if i > 1: mz.migration.gas[i - 1][i] = 1e-3 * (i + 1)
		if i > 61: 
			mz.zones[i].tau_star = float("inf") 
			for j in mz.zones[i].elements: 
				mz.zones[i].entrainment.agb[j] = 0 
				mz.zones[i].entrainment.ccsne[j] = 0 
				mz.zones[i].entrainment.sneia[j] = 0 
		# else: 
		# 	mz.zones[i].tau_star = 2 * ((RAD_BINS[i + 1]**2 - RAD_BINS[i]**2) / 
		# 		RAD_BINS[1]**2 * m.exp((-RAD_BINS[i] + 0.125) / 
		# 			float(sys.argv[1])))**(-0.5) 
	print("Running....") 
	mz.run(np.linspace(0, 13.8, 1381), overwrite = True) 

if __name__ == "__main__": 
	run_simulation() 

