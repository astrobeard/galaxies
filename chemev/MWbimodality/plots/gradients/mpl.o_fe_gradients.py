""" 
Plots the [O/H] and [Fe/H] gradients of both gas and tracer particles for a 
given multizone output. 

ARGV 
==== 
1) 		The name of the multizone output 
2) 		The name of the output image 
""" 

import matplotlib.pyplot as plt 
import plots 
plots.mpltoolkit.load_mpl_presets() 
from scipy import stats
import numpy as np 
import math as m 
import vice 
import sys 
import os 

def setup_axes(): 
	fig = plt.figure(figsize = (14, 7)) 
	axes = 2 * [None] 
	for i in range(2): 
		axes[i] = fig.add_subplot(121 + i, facecolor = "white") 
		axes[i].set_xlabel(r"$R_\text{gal}$ [kpc]") 
		axes[i].set_xlim([-2, 32]) 
	axes[0].set_ylim([-1.1, 0.7]) 
	axes[0].yaxis.set_ticks([-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6])  
	axes[1].set_ylim([-0.07, 0.27]) 
	axes[0].set_ylabel("[X/H]") 
	axes[1].set_ylabel("[O/Fe]") 
	return axes 


# def plot_tracers(axes, multiout): 
# 	tracers = [list(i) for i in zip(
# 		multiout.tracers["zone_final"], 
# 		multiout.tracers["mass"], 
# 		multiout.tracers["z(o)"], 
# 		multiout.tracers["z(fe)"], 
# 	)] 
# 	O = len(multiout.zones.keys()) * [0.] 
# 	Fe = len(multiout.zones.keys()) * [0.] 
# 	OFe = len(multiout.zones.keys()) * [0.] 
# 	for i in range(len(multiout.zones.keys())): 
# 		rgal = 0.25 * i + 0.125 
# 		inzone = list(filter(lambda x: x[0] == i, tracers)) 
# 		# try: 
# 		# 	O[i] = m.log10(stats.mode([row[2] for row in inzone])[0][0] / 
# 		# 		vice.solar_z["o"]) 
# 		# except ValueError: 
# 		# 	O[i] = -float("inf") 
# 		# try: 
# 		# 	Fe[i] = m.log10(stats.mode([row[3] for row in inzone])[0][0] / 
# 		# 		vice.solar_z["fe"]) 
# 		# except ValueError: 
# 		# 	Fe[i] = -float("inf") 
# 		inzone = list(filter(lambda x: x[2] != 0 and x[3] != 0, inzone)) 
# 		O[i] = m.log10(np.median([row[2] for row in inzone]) / 
# 			vice.solar_z["o"]) 
# 		Fe[i] = m.log10(np.median([row[3] for row in inzone]) / 
# 			vice.solar_z["fe"]) 
# 		OFe_ = [m.log10(i[2] / vice.solar_z["o"]) - m.log10(i[3] / 
# 			vice.solar_z["fe"]) for i in inzone] 
# 		OFe[i] = np.median(OFe_) 

# 	radii = [0.25 * i for i in range(len(multiout.zones.keys()))] 
# 	axes[0].scatter(radii, O, c = plots.mpltoolkit.named_colors()["red"], 
# 		marker = plots.mpltoolkit.markers()["square"]) 
# 	axes[0].scatter(radii, Fe, c = plots.mpltoolkit.named_colors()["blue"], 
# 		marker = plots.mpltoolkit.markers()["square"]) 
# 	axes[1].scatter(radii, OFe, c = plots.mpltoolkit.named_colors()["black"], 
# 		marker = plots.mpltoolkit.markers()["square"]) 

def mode_stellar_metallicity(zone, mdf_key): 
	try: 
		idx = zone.mdf[mdf_key].index(max(zone.mdf[mdf_key])) 
		return (zone.mdf["bin_edge_left"][idx] + 
			zone.mdf["bin_edge_right"][idx]) / 2. 
	except ValueError: 
		return float("nan") 


def stellar_dispersion(zone, mdf_key): 
	s = 0 
	low = 0 
	high = 0 
	for i in range(len(zone.mdf[mdf_key])): 
		s += (zone.mdf["bin_edge_right"][i] - 
			zone.mdf["bin_edge_left"][i]) * zone.mdf[mdf_key][i]
		if s >= 0.159 and low == 0: 
			low = (zone.mdf["bin_edge_left"][i] + 
				zone.mdf["bin_edge_right"][i]) / 2. 
		if s >= 0.841: 
			high = (zone.mdf["bin_edge_left"][i] + 
				zone.mdf["bin_edge_right"][i]) / 2. 
			break 
	return [low, high] 


def plot_stellar_metallicities(axes, multiout): 
	O = len(multiout.zones.keys()) * [0.] 
	O_disp = len(multiout.zones.keys()) * [0.] 
	Fe = len(multiout.zones.keys()) * [0.] 
	Fe_disp = len(multiout.zones.keys()) * [0.] 
	OFe = len(multiout.zones.keys()) * [0.] 
	OFe_disp = len(multiout.zones.keys()) * [0.] 
	for i in range(len(multiout.zones.keys())): 
		O[i] = mode_stellar_metallicity(multiout.zones["zone%d" % (i)], 
			"dn/d[o/h]") 
		O_disp[i] = stellar_dispersion(multiout.zones["zone%d" % (i)], 
			"dn/d[o/h]") 
		Fe[i] = mode_stellar_metallicity(multiout.zones["zone%d" % (i)], 
			"dn/d[fe/h]") 
		Fe_disp[i] = stellar_dispersion(multiout.zones["zone%d" % (i)], 
			"dn/d[fe/h]") 
		OFe[i] = mode_stellar_metallicity(multiout.zones["zone%d" % (i)], 
			"dn/d[o/fe]") 
		OFe_disp[i] = stellar_dispersion(multiout.zones["zone%d" % (i)], 
			"dn/d[o/fe]") 
	radii = [0.25 * i for i in range(len(multiout.zones.keys()))] 
	axes[0].fill_between(radii, [row[0] for row in O_disp], 
		[row[1] for row in O_disp], 
		color = plots.mpltoolkit.named_colors()["red"], 
		alpha = 0.3) 
	axes[0].fill_between(radii, [row[0] for row in Fe_disp], 
		[row[1] for row in Fe_disp], 
		color = plots.mpltoolkit.named_colors()["blue"], 
		alpha = 0.3)  
	axes[1].fill_between(radii, [row[0] for row in OFe_disp], 
		[row[1] for row in OFe_disp], 
		color = plots.mpltoolkit.named_colors()["black"], 
		alpha = 0.3) 

	axes[0].scatter(radii, O, c = plots.mpltoolkit.named_colors()["red"], 
		marker = plots.mpltoolkit.markers()["star"], s = 50) 
	axes[0].scatter(radii, Fe, c = plots.mpltoolkit.named_colors()["blue"], 
		marker = plots.mpltoolkit.markers()["star"], s = 50) 
	axes[1].scatter(radii, OFe, c = plots.mpltoolkit.named_colors()["black"], 
		marker = plots.mpltoolkit.markers()["star"], s = 50) 

def plot_gas_phase_metallicities(axes, multiout): 
	O = [multiout.zones["zone%d" % (i)].history["[o/h]"][-1] for i in range(
		len(multiout.zones.keys()))] 
	Fe = [multiout.zones["zone%d" % (i)].history["[fe/h]"][-1] for i in range(
		len(multiout.zones.keys()))] 
	OFe = [multiout.zones["zone%d" % (i)].history["[o/fe]"][-1] for i in range(
		len(multiout.zones.keys()))] 
	radii = [0.25 * i for i in range(len(multiout.zones.keys()))] 
	axes[0].plot(radii, O, c = plots.mpltoolkit.named_colors()["red"]) 
		# marker = plots.mpltoolkit.markers()["circle"]) 
	axes[0].plot(radii, Fe, c = plots.mpltoolkit.named_colors()["blue"]) 
		# marker = plots.mpltoolkit.markers()["circle"]) 
	axes[1].plot(radii, OFe, c = plots.mpltoolkit.named_colors()["black"]) 
		# marker = plots.mpltoolkit.markers()["circle"]) 


def O_Fe_legend(ax): 
	O = ax.scatter(0, 0, c = plots.mpltoolkit.named_colors()["white"], 
		label = "O") 
	Fe = ax.scatter(0, 0, c = plots.mpltoolkit.named_colors()["white"], 
		label = "Fe") 
	leg = ax.legend(loc = plots.mpltoolkit.mpl_loc("upper right"), 
		ncol = 1, frameon = False, bbox_to_anchor = (0.98, 0.98)) 
	O.remove() 
	Fe.remove() 
	for i in range(2): 
		leg.get_texts()[i].set_color(["red", "blue"][i]) 


if __name__ == "__main__": 
	plt.clf() 
	axes = setup_axes() 
	out = vice.multioutput(sys.argv[1]) 
	plot_stellar_metallicities(axes, out) 
	plot_gas_phase_metallicities(axes, out) 
	O_Fe_legend(axes[0]) 
	plt.tight_layout() 
	plt.savefig(sys.argv[2]) 
	plt.clf() 





