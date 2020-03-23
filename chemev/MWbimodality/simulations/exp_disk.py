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

TIME_BINS = np.linspace(0, 12.8, 41).tolist() 
RAD_BINS = np.linspace(0, 30, 121).tolist() 
ZONE_WIDTH = 0.25 


def eta(rgal): 
	""" 
	The relation between the mass loading factor eta and galactocentric radius 

	rgal :: galactocentric radius in kpc 

	Notes 
	===== 
	0.00185 is the Mg CCSN yield assuming our O yield and [O/Mg] = 0 
	-0.6 is r - 1 where r is the instantaneous recycling parameter from WAF17 
	""" 
	return 0.00185 / vice.solar_z["mg"] * 10**(0.06 * (rgal - 4) - 0.3) - 0.6 


def run_simulation(): 
	from vice.yields.presets import my_yields 
	mz = vice.multizone(name = "moddisk_vigorousSF", 
		n_zones = len(RAD_BINS) - 1, 
		n_tracers = 2, verbose = True, simple = False) 
	mz.migration.stars = tracers.UWhydro(TIME_BINS, RAD_BINS, 
		filename = "%s_extra_tracer_data.out" % (mz.name)) 
	for i in range(mz.n_zones): 
		mz.zones[i].mode = "gas" 
		mz.zones[i].func = gas_disks.static_exponential(i, 6.0e9, 
			RAD_BINS, float(sys.argv[1]))  
		mz.zones[i].bins = np.linspace(-3, 1, 401) 
		mz.zones[i].elements = ["mg", "fe", "o"] 
		mz.zones[i].eta = eta( (RAD_BINS[i] + RAD_BINS[i + 1]) / 2 ) 
		mz.zones[i].dt = 0.02 
		if i > 61: 
			mz.zones[i].tau_star = float("inf") 
			for j in mz.zones[i].elements: 
				mz.zones[i].entrainment.agb[j] = 0 
				mz.zones[i].entrainment.ccsne[j] = 0 
				mz.zones[i].entrainment.sneia[j] = 0 
		else: 
			# tau_star ~ e^r/(2r_s) ; extra factor of two for bin center 
			mz.zones[i].tau_star = 0.2 * m.exp( (RAD_BINS[i] + RAD_BINS[i + 1]) / 
				(4 * float(sys.argv[1]))) 
	print("Running....") 
	mz.run(np.linspace(0, 12.8, 641), overwrite = True) 
	mz.migration.stars.close_file() 

if __name__ == "__main__": 
	run_simulation() 

