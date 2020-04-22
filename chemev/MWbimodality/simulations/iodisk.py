r""" 
Runs a simulation of an inside-out disk with VICE's multizone object with 
tracer particles tuned to the UWhydro simulation. 

ARGV 
----
1) 		The name of the output 
""" 

import tracers 
import gas_disks 
from common import * 
import numpy as np 
import math as m 
import vice 
from vice.yields.presets import my_yields 
import sys 
import os 

TAU_STAR0 = 0.2

def tau_in(rgal): 
	r""" 
	The mathematical relation between galactocentric radius and the infall 
	timescale. 

	Parameters 
	----------
	rgal : real number 
		Galactocentric radius in kpc. 

	Returns 
	-------
	tau_in : real number 
		The infall timescale in Gyr. 
	""" 
	return 1 + (rgal + 1.e-12) / 1.5 
	# return 6 


def Min0(rgal, k = 0.1, tau_star0 = TAU_STAR0, scale = 3): 
	t_star = tau_star(rgal, norm = tau_star0, scale = scale) 
	t_in = tau_in(rgal) 
	return k * (
		rgal / ZONE_WIDTH * m.exp(-rgal / scale) * 
		harmonic_timescale(t_in, depletion_time(t_star, 
			eta = eta(rgal, corrective = t_star / t_in))
			# eta = eta(rgal)) 
		)**(-1) * (
			m.exp(-12.8 / depletion_time(t_star, 
				eta = eta(rgal, corrective = t_star / t_in))) - 
				# eta = eta(rgal))) - 
			m.exp(-12.8 / t_in) 
		)**(-1) 
	)


def sfr_norm(r, rs = 3, k = 100): 
	return k * tau_in(r)**(-2) * (1 - (1 + 12.8 / tau_in(r)) * m.exp(-12.8 / 
		tau_in(r)))**(-1) * 2 * m.pi * r * m.exp(-r / rs) * 0.25 


def lintexp_sfr_norm(r, rs = 3, k = 1000): 
	t = 12.8 
	t1 = 1 
	return k * r * m.exp(-r / rs) / t * (0.5 * t + tau_in(r) * (1 - 
		m.exp(-(t - t1) / tau_in(r))))**(-1) 


def run_simulation(): 
	mz = vice.multizone(name = sys.argv[1], n_zones = len(RAD_BINS) - 1, 
		n_stars = 4, verbose = True, simple = False) 
	mz.migration.stars = tracers.UWhydro(TIME_BINS, RAD_BINS, 
		n_stars = mz.n_stars, 
		filename = "%s_extra_tracer_data.out" % (mz.name)) 
	# mz.migration.stars = tracers.UWhydro_1event(TIME_BINS, RAD_BINS, 
	# 	n_stars = mz.n_stars, 
		# filename = "%s_extra_tracer_data.out" % (mz.name)) 
	for i in range(mz.n_zones): 
		# mz.zones[i].mode = "ifr" 
		mz.zones[i].mode = "sfr" 
		# mz.zones[i].func = gas_disks.exponential_decay(
		# 	Min0( (RAD_BINS[i] + RAD_BINS[i + 1]) / 2 ), 
		# 	tau_in(RAD_BINS[i] + RAD_BINS[i + 1]) / 2) 
		# mz.zones[i].func = gas_disks.linear_exponential(
		#  	sfr_norm((RAD_BINS[i] + RAD_BINS[i + 1]) / 2), 
		#  	tau_in((RAD_BINS[i] + RAD_BINS[i + 1]) / 2))  
		mz.zones[i].func = gas_disks.linear_then_exponential(
			lintexp_sfr_norm((RAD_BINS[i] + RAD_BINS[i + 1]) / 2), 
			tau_in((RAD_BINS[i] + RAD_BINS[i + 1]) / 2), 
			1) 
		mz.bins = np.linspace(-3, 1, 401) 
		mz.zones[i].elements = ["mg", "fe", "o"] 
		mz.zones[i].dt = 0.01 
		mz.zones[i].Mg0 = 0 
		if i > 61: 
			mz.zones[i].func = lambda t: 0 
			mz.zones[i].tau_star = 100 
			# mz.zones[i].tau_star = float("inf") 
			mz.zones[i].eta = 100 
			for j in mz.zones[i].elements: 
				mz.zones[i].entrainment.agb[j] = 0 
				mz.zones[i].entrainment.ccsne[j] = 0 
				mz.zones[i].entrainment.sneia[j] = 0 
		else: 
			mz.zones[i].tau_star = tau_star(ZONE_WIDTH * (i + 0.5), 
				norm = TAU_STAR0) 
			mz.zones[i].eta = eta( (RAD_BINS[i] + RAD_BINS[i + 1]) / 2 , 
				corrective = mz.zones[i].tau_star / tau_in(
					(RAD_BINS[i] + RAD_BINS[i + 1]) / 2)) 
			# mz.zones[i].eta = eta( (RAD_BINS[i] + RAD_BINS[i + 1]) / 2 ) 
		mz.zones[i].schmidt = True 
	print("Running....") 
	# mz.run(np.linspace(0, 12.8, 641), overwrite = True) 
	mz.run(np.linspace(0, 12.8, 257), overwrite = True) 
	mz.migration.stars.close_file() 


if __name__ == "__main__": 
	# for i in RAD_BINS[:60]: 
	# 	print("R = %.2f kpc ; Min0 = %.2f ; tau_in = %.2f" % (
	# 		i + ZONE_WIDTH / 2, Min0(i + ZONE_WIDTH / 2), 
	# 		tau_in(i + ZONE_WIDTH / 2))) 
	# for i in RAD_BINS[:60]: 
	# 	print("R = %.2f kpc ; sfr_norm = %.2f ; tau_sfh = %.2f" % (
	# 		i + ZONE_WIDTH / 2, sfr_norm(i + ZONE_WIDTH / 2), 
	# 		tau_in(i + ZONE_WIDTH / 2))) 
	for i in RAD_BINS[:60]: 
		print("R = %.2f kpc; norm = %.2f ; tau_sfh = %.2f" % (
			i + ZONE_WIDTH / 2, lintexp_sfr_norm(i + ZONE_WIDTH / 2), 
			tau_in(i + ZONE_WIDTH / 2))) 
	run_simulation() 

