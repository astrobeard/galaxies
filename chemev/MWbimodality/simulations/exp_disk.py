""" 
Runs a simulation of an exponential disk with VICE's multizone object 
assuming constant SFR with tracer particles tuned to given output data 

ARGV 
==== 
1) 	The scale length of the disk in kpc 
""" 

import tracers 
import gas_disks 
from common import * 
import numpy as np 
import math as m 
import vice 
import sys 
import os 


def run_simulation(): 
	from vice.yields.presets import my_yields 
	mz = vice.multizone(name = "moddisk", 
		n_zones = len(RAD_BINS) - 1, 
		n_tracers = 4, verbose = True, simple = False) 
	mz.migration.stars = tracers.UWhydro(TIME_BINS, RAD_BINS, 
		filename = "%s_extra_tracer_data.out" % (mz.name)) 
	for i in range(mz.n_zones): 
		mz.zones[i].mode = "gas" 
		mz.zones[i].func = gas_disks.static_exponential(i, 6.0e9, 
			RAD_BINS, float(sys.argv[1]))  
		mz.zones[i].bins = np.linspace(-3, 1, 401) 
		mz.zones[i].elements = ["mg", "fe", "o"] 
		mz.zones[i].eta = eta( (RAD_BINS[i] + RAD_BINS[i + 1]) / 2 ) 
		mz.zones[i].dt = 0.01 
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
	print("Running....") 
	mz.run(np.linspace(0, 12.8, 641), overwrite = True) 
	mz.migration.stars.close_file() 

if __name__ == "__main__": 
	run_simulation() 

