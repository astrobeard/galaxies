""" 
Plots the age-metallicity relation for stars and gas color-coded by radius (of 
birth in the case of stars) for stars in a given annulus. 

ARGV 
==== 
1) 		The name of the VICE output 
2) 		The name of the output image
""" 

import matplotlib.pyplot as plt 
import plots 
plots.mpltoolkit.load_mpl_presets() 
import numpy as np 
import math as m 
import vice 
import sys 
import os 

CMAP = "plasma" 
zone_min = int(7 / 0.25) 
zone_max = int(9 / 0.25) 

def setup_axes(): 
	fig = plt.figure(figsize = (21, 7)) 
	axes = 3 * [None] 
	for i in range(3): 
		axes[i] = fig.add_subplot(131 + i, facecolor = "white") 
		axes[i].set_xlabel("Age [Gyr]") 
	axes[0].set_ylabel("[O/H]") 
	axes[1].set_ylabel("[Fe/H]") 
	axes[2].set_ylabel("[O/Fe]") 
	return axes 

def plot_tracers(axes, multiout): 
	cmap = plt.get_cmap(CMAP) 
	tracers = [list(i) for i in zip(
		multiout.tracers["formation_time"], 
		multiout.tracers["zone_final"], 
		multiout.tracers["zone_origin"], 
		multiout.tracers["mass"], 
		multiout.tracers["z(O)"], 
		multiout.tracers["z(fe)"] 
	)]
	tracers = list(filter(lambda x: zone_min <= x[1] <= zone_max, 
		tracers)) 
	ages = len(tracers) * [0.] 
	sizes = len(tracers) * [0.] 
	colors = len(tracers) * [0.] 
	OH = len(tracers) * [0.] 
	FeH = len(tracers) * [0.] 
	OFe = len(tracers) * [0.] 
	for i in range(len(tracers)): 
		ages[i] = 13.8 - tracers[i][0] 
		OH[i] = m.log10(tracers[i][4] / vice.solar_z['o']) 
		FeH[i] = m.log10(tracers[i][5] / vice.solar_z['fe']) 
		OFe[i] = OH[i] - FeH[i] 
		sizes[i] = tracers[i][3] / 4e6 * 4 * (1 - 
			vice.cumulative_return_fraction(tracers[i][0])) 
		colors[i] = 0.25 * tracers[i][2] 
	sc = axes[0].scatter(ages, OH, c = colors, s = sizes, cmap = cmap, 
		vmin = 0, vmax = 15) 
	return sc 


if __name__ == "__main__": 
	plt.clf() 
	axes = setup_axes() 
	out = vice.multioutput(sys.argv[1]) 
	plot_tracers(axes, out) 
	plt.tight_layout() 
	plt.savefig(sys.argv[2]) 
	plt.clf() 



