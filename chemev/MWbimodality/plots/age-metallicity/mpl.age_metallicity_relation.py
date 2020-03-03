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

CMAP = "plasma_r" 
zone_min = int(7 / 0.25) 
zone_max = int(8.75 / 0.25) 

def setup_axes(): 
	fig = plt.figure(figsize = (21, 7)) 
	axes = 3 * [None] 
	for i in range(3): 
		axes[i] = fig.add_subplot(131 + i, facecolor = "white") 
		axes[i].set_xlabel("Age [Gyr]") 
		axes[i].set_xlim([-1, 15]) 
	axes[0].set_ylabel("[O/H]") 
	axes[1].set_ylabel("[Fe/H]") 
	axes[2].set_ylabel("[O/Fe]") 
	axes[0].set_ylim([-0.7, 0.5]) 
	axes[1].set_ylim([-0.7, 0.5]) 
	axes[2].set_ylim([0.0, 0.5]) 
	return axes 

def plot_tracers(axes, tracers): 
	cmap = plt.get_cmap(CMAP) 
	tracers = [list(i) for i in zip(
		tracers["formation_time"], 
		tracers["zone_final"], 
		tracers["zone_origin"], 
		tracers["mass"], 
		tracers["z(o)"], 
		tracers["z(fe)"] 
	)]
	ages = len(tracers) * [0.] 
	sizes = len(tracers) * [0.] 
	colors = len(tracers) * [0.] 
	OH = len(tracers) * [0.] 
	FeH = len(tracers) * [0.] 
	OFe = len(tracers) * [0.] 
	for i in range(len(tracers)): 
		ages[i] = 13.8 - tracers[i][0] 
		if tracers[i][4]: 
			OH[i] = m.log10(tracers[i][4] / vice.solar_z['o']) 
		else: 
			OH[i] = -float("inf") 
		if tracers[i][5]: 
			FeH[i] = m.log10(tracers[i][5] / vice.solar_z['fe']) 
		else: 
			FeH[i] = -float("inf") 
		OFe[i] = OH[i] - FeH[i] 
		sizes[i] = tracers[i][3] / 4e7 * 4 * (1 - 
			vice.cumulative_return_fraction(tracers[i][0])) 
		colors[i] = 0.25 * tracers[i][2] 
	axes[0].scatter(ages, OH, c = colors, s = sizes, cmap = cmap, 
		vmin = 0, vmax = 15) 
	axes[1].scatter(ages, FeH, c = colors, s = sizes, cmap = cmap, 
		vmin = 0, vmax = 15) 
	sc = axes[2].scatter(ages, OFe, c = colors, s = sizes, cmap = cmap, 
		vmin = 0, vmax = 15) 
	return sc 


def weighted_median(values, weights, stop = 0.5): 
	indeces = np.argsort(values) 
	values = [values[i] for i in indeces] 
	weights = [weights[i] for i in indeces] 
	weights = [i / sum(weights) for i in weights] 
	s = 0 
	for i in range(len(weights)): 
		s += weights[i] 
		if s > stop: 
			idx = i - 1 
			break 
	return values[idx] 


def feuillet_points(ax, tracers, element): 
	bins = np.arange(-1, 1.1, 0.1) 
	ages = (len(bins) - 1) * [0.] 
	lowers = (len(bins) - 1) * [0.] 
	uppers = (len(bins) - 1) * [0.] 
	for i in range(len(ages)): 
		stars = tracers.filter("[%s/h]" % (element), ">=", bins[i]) 
		stars = stars.filter("[%s/h]" % (element), "<=", bins[i + 1]) 
		if len(stars["age"]) > 20: 
			""" 
			Calculate the mass-weighted median 
			""" 
			masses = list(map(lambda x, y: x * (1 - 
				vice.cumulative_return_fraction(y)), 
				stars["mass"], stars["age"])) 
			ages[i] = weighted_median(stars["age"], masses) 
			lowers[i] = weighted_median(stars["age"], masses, stop = 0.16) 
			uppers[i] = weighted_median(stars["age"], masses, stop = 0.84) 
		else: 
			ages[i] = float("nan") 
			lowers[i] = float("nan") 
			uppers[i] = float("nan") 
	ax.scatter(ages, list(map(lambda x, y: (x + y) / 2, bins[1:], bins[:-1])), 
		marker = plots.mpltoolkit.markers()["star"], 
		c = plots.mpltoolkit.named_colors()["black"], s = 100) 
	ax.errorbar(list(map(lambda x, y: (x + y) / 2, uppers, lowers)), 
		list(map(lambda x, y: (x + y) / 2, bins[1:], bins[:-1])), 
		xerr = list(map(lambda x, y: (x - y) / 2, uppers, lowers)), 
		c = plots.mpltoolkit.named_colors()["black"], 
		linestyle = "None")  


if __name__ == "__main__": 
	plt.clf() 
	axes = setup_axes() 
	out = vice.multioutput(sys.argv[1]) 
	extra_tracer_data = np.genfromtxt("%s_extra_tracer_data.out" % (out.name)) 
	out.tracers["zfinal"] = [row[-1] for row in extra_tracer_data[:out.tracers.size[0]]] 
	fltrd_tracers = out.tracers.filter("zfinal", ">=", -3.) 
	fltrd_tracers = fltrd_tracers.filter("zfinal", "<=", 3.) 
	fltrd_tracers = fltrd_tracers.filter("zone_final", ">=", zone_min) 
	fltrd_tracers = fltrd_tracers.filter("zone_final", "<=", zone_max) 
	sc = plot_tracers(axes, fltrd_tracers) 
	cbar = plt.colorbar(sc, 
		cax = plots.mpltoolkit.append_axes(axes[2]), pad = 0.0) 
	cbar.set_label(r"$R_\text{gal}$ of birth [kpc]") 
	feuillet_points(axes[0], fltrd_tracers, "o") 
	feuillet_points(axes[1], fltrd_tracers, "fe") 
	plt.tight_layout() 
	plt.savefig(sys.argv[2]) 
	plt.clf() 



