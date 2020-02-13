""" 
Runs a simulation of an exponential disk with VICE's multizone object 
assuming constant SFR with tracer particles tuned to given output data 

ARGV 
==== 
1) 	The scale length of the disk in kpc 
""" 

import tracers 
import gas_disks 
import numpy as np 
import math as m 
import vice 
import sys 
import os 

TIME_BINS = np.linspace(0, 13.8, 21).tolist() 
RAD_BINS = np.linspace(0, 30, 121).tolist() 
ZONE_WIDTH = 0.25 


def run_simulation(): 
	from vice.yields.presets import my_yields 
	mz = vice.multizone(name = "expdisk_schmidt_zfilter_lowtaustar", 
		n_zones = len(RAD_BINS) - 1, 
		n_tracers = 2, verbose = True, simple = False) 
	mz.migration.stars = tracers.UWhydro_zfilter(TIME_BINS, RAD_BINS) 
	for i in range(mz.n_zones): 
		mz.zones[i].mode = "gas" 
		mz.zones[i].func = gas_disks.static_exponential(i, 6.0e9, 
			RAD_BINS, float(sys.argv[1]))  
		mz.zones[i].bins = np.linspace(-3, 1, 401) 
		mz.zones[i].elements = ["mg", "fe", "o"] 
		mz.zones[i].eta = 0.2 * (i + 1) 
		mz.zones[i].dt = 0.02 
		# mz.zones[i].schmidt = True 
		# if i > 1: mz.migration.gas[i - 1][i] = 1e-3 * (i + 1)
		if i > 61: 
			mz.zones[i].tau_star = float("inf") 
			for j in mz.zones[i].elements: 
				mz.zones[i].entrainment.agb[j] = 0 
				mz.zones[i].entrainment.ccsne[j] = 0 
				mz.zones[i].entrainment.sneia[j] = 0 
		else: 
			# tau_star ~ e^r/(2r_s) ; extra factor of two for bin center 
			mz.zones[i].tau_star = 2 * m.exp( (RAD_BINS[i] + RAD_BINS[i + 1]) / 
				(4 * float(sys.argv[1]))) 
		# else: 
		# 	mz.zones[i].tau_star = 2 * ((RAD_BINS[i + 1]**2 - RAD_BINS[i]**2) / 
		# 		RAD_BINS[1]**2 * m.exp((-RAD_BINS[i] + 0.125) / 
		# 			float(sys.argv[1])))**(-0.5) 
	print("Running....") 
	mz.run(np.linspace(0, 13.8, 691), overwrite = True) 

if __name__ == "__main__": 
	run_simulation() 

