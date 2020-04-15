""" 
Produces a plot of the metallicity gradient for a given multizone output 

ARGV 
----
1) 		The name of the multizone output 
2) 		The name of the output image (no extension) 
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
	""" 
	Set up the matplotlib axis to plot on 
	""" 
	fig = plt.figure(figsize = (14, 7)) 
	axes = 2 * [None] 
	for i in range(2): 
		axes[i] = fig.add_subplot(121 + i, facecolor = "white") 
		axes[i].set_xlabel(r"$R_\text{gal}$ [kpc]") 
		axes[i].set_xlim([-1, 21]) 
	axes[0].set_ylim([-0.55, 0.75]) 
	axes[1].set_ylim([-0.1, 0.5]) 
	axes[0].set_ylabel("[X/H]") 
	axes[1].set_ylabel("[O/Fe]") 
	return axes 


def mode_stellar_metallicity(zone, mdf_key): 
	""" 
	Find the mode of the stellar metallicity distribution 

	zone : vice.output 
		The output object from an individual zone in the multizone simulation 
	mdf_key : str 
		The string for the mdf dataframe 
	""" 
	try: 
		idx = zone.mdf[mdf_key].index(max(zone.mdf[mdf_key])) 
		return (zone.mdf["bin_edge_left"][idx] + 
			zone.mdf["bin_edge_right"][idx]) / 2. 
	except ValueError: 
		return float("nan") 


def stellar_dispersion(zone, mdf_key): 
	""" 
	Find the 16th and 84th percentile of the metallicity distribution of 
	stars 

	zone : vice.output 
		The output object from an individual zone in the multizone simulation 
	mdf_key : str 
		The string for the mdf dataframe 
	""" 
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


def plot_stellar_metallicities(axes, multioutput): 
	""" 
	Plot the stellar metallicity information 

	axes : list 
		A 2-element list of matplotlib subplots to put [X/H] v. radius on the 
		first and [O/Fe] v. radius on the second 
	multioutput : vice.multioutput 
		The multioutput object from the simulation 
	""" 
	O = [mode_stellar_metallicity(multioutput.zones["zone%d" % (i)], 
		"dn/d[o/h]") for i in range(len(multioutput.zones.keys()))] 
	O_disp = [stellar_dispersion(multioutput.zones["zone%d" % (i)], 
		"dn/d[o/h]") for i in range(len(multioutput.zones.keys()))] 
	Fe = [mode_stellar_metallicity(multioutput.zones["zone%d" % (i)], 
		"dn/d[fe/h]") for i in range(len(multioutput.zones.keys()))] 
	Fe_disp = [stellar_dispersion(multioutput.zones["zone%d" % (i)], 
		"dn/d[fe/h]") for i in range(len(multioutput.zones.keys()))] 
	OFe = [mode_stellar_metallicity(multioutput.zones["zone%d" % (i)], 
		"dn/d[o/fe]") for i in range(len(multioutput.zones.keys()))] 
	OFe_disp = [stellar_dispersion(multioutput.zones["zone%d" % (i)], 
		"dn/d[o/fe]") for i in range(len(multioutput.zones.keys()))] 
	radii = [0.25 * i + 0.125 for i in range(len(multioutput.zones.keys()))] 
	axes[0].scatter(radii, O, c = plots.mpltoolkit.named_colors()["red"], 
		marker = plots.mpltoolkit.markers()["star"], s = 50, zorder = 20) 
	axes[0].scatter(radii, Fe, c = plots.mpltoolkit.named_colors()["blue"], 
		marker = plots.mpltoolkit.markers()["star"], s = 50, zorder = 20) 
	axes[1].scatter(radii, OFe, c = plots.mpltoolkit.named_colors()["black"], 
		marker = plots.mpltoolkit.markers()["star"], s = 50, zorder = 20) 
	axes[0].fill_between(radii, [row[0] for row in O_disp], 
		[row[1] for row in O_disp], alpha = 0.3, zorder = 0, 
		color = plots.mpltoolkit.named_colors()["red"])
	axes[0].fill_between(radii, [row[0] for row in Fe_disp], 
		[row[1] for row in Fe_disp], alpha = 0.3, zorder = 0, 
		color = plots.mpltoolkit.named_colors()["blue"]) 
	axes[1].fill_between(radii, [row[0] for row in OFe_disp], 
		[row[1] for row in OFe_disp], alpha = 0.3, zorder = 0, 
		color = plots.mpltoolkit.named_colors()["black"]) 


def plot_gas_phase_metallicity(axes, out): 
	""" 
	Plot the gas-phase metallicities as a function of galactocentric radius 

	axes : list 
		A 2-element list of matplotlib subplots to put [X/H] v. radius on 
		the first and [O/Fe] v. radius on the second 
	out : vice.multioutput 
		The multioutput object from the simulation  
	""" 
	O = [out.zones["zone%d" % (i)].history["[o/h]"][-1] for i in range(
		len(out.zones.keys()))] 
	Fe = [out.zones["zone%d" % (i)].history["[fe/h]"][-1] for i in range(
		len(out.zones.keys()))] 
	OFe = [out.zones["zone%d" % (i)].history["[o/fe]"][-1] for i in range(
		len(out.zones.keys()))] 
	radii = [0.25 * i + 0.125 for i in range(len(out.zones.keys()))] 
	axes[0].plot(radii[:62], O[:62], 
		c = plots.mpltoolkit.named_colors()["red"]) 
	axes[0].plot(radii[:62], Fe[:62], 
		c = plots.mpltoolkit.named_colors()["blue"]) 
	axes[1].plot(radii[:62], OFe[:62], 
		c = plots.mpltoolkit.named_colors()["black"]) 


def legend(ax): 
	""" 
	Put a red 'O' and blue 'Fe' in the upper right of the subplot ax 
	""" 
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
	output = vice.multioutput(sys.argv[1]) 
	plot_stellar_metallicities(axes, output) 
	plot_gas_phase_metallicity(axes, output) 
	legend(axes[0]) 
	plt.tight_layout() 
	plt.savefig("%s.pdf" % (sys.argv[2])) 
	plt.savefig("%s.png" % (sys.argv[2])) 
	plt.clf() 

