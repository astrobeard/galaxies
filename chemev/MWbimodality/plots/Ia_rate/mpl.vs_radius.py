""" 
This script produces a plot of the Ia rate in each annular zone in comparison to what 
is expected given the star formation history 

ARGV 
==== 
1) 	The name of the multizone output 
2)	The name of the output image 
""" 

import matplotlib.pyplot as plt 
import plots 
plots.mpltoolkit.load_mpl_presets() 
import numpy as np 
import math as m 
import vice 
import sys 
import os 

XLIM = [8, 18] 

def setup_axes(): 
	fig = plt.figure(figsize = (7, 7)) 
	ax = fig.add_subplot(111, facecolor = "white") 
	ax.set_ylabel(r"$\dot{N}_\text{Ia}$ [yr$^{-1}$] (unnormalized)") 
	bottom = plots.mpltoolkit.append_axes(ax, size = "20%", loc = "bottom", 
		sharex = ax) 
	ax.set_xlim(XLIM) 
	ax.set_ylim([-0.02, 0.22]) 
	bottom.set_ylim([-0.01, 0.03]) 
	plt.setp(ax.get_xticklabels(), visible = False) 
	bottom.xaxis.set_ticks([8, 10, 12, 14, 16, 18]) 
	bottom.yaxis.set_ticks([0.0, 0.02])
	bottom.set_xlabel(r"$R_\text{gal}$ [kpc]") 
	bottom.set_ylabel(r"$\Delta$") 
	return [ax, bottom] 


def R_Ia(t): 
	if t < 0.15: 
		return 0 
	else: 
		return t**-1.1 


def get_expected(out): 
	rates = len(out.zones.keys()) * [0.] 
	for i in range(len(rates)): 
		for j in range(len(out.zones["zone%d" % (i)].history["time"])): 
			rates[i] += R_Ia(
				13.8 - out.zones["zone%d" % (i)].history["time"][j] 
			) * out.zones["zone%d" % (i)].history["sfr"][j] * 1e9 * (
				out.zones["zone%d" % (i)].history["time"][1] - 
				out.zones["zone%d" % (i)].history["time"][0] 
			) 
		rates[i] /= 1e11 
	return rates 
	# radii = [0.25 * i + 0.125 for i in range(len(out.zones.keys()))] 
	# ax.plot(radii, rates, c = plots.mpltoolkit.named_colors()["black"], 
	# 	linestyle = '--') 


def get_actual(out): 
	rates = len(out.zones.keys()) * [0.] 
	for i in range(len(out.zones.keys())): 
		stars = out.tracers.filter("zone_final", "==", i) 
		for j in range(len(stars["mass"])): 
			rates[i] += stars["mass"][j] * R_Ia(
				13.8 - stars["formation_time"][j]
			) 
		rates[i] /= 1e11 
	return rates 
	# radii = [0.25 * i + 0.125 for i in range(len(out.zones.keys()))] 
	# ax.scatter(radii, rates, marker = plots.mpltoolkit.markers()["x"], 
	# 	color = plots.mpltoolkit.named_colors()["black"]) 


def plot_rates(ax, bottom, out): 
	radii = [0.25 * i + 0.125 for i in range(len(out.zones.keys()))] 
	expected = get_expected(out) 
	actual = get_actual(out) 
	ax.plot(radii, expected, c = plots.mpltoolkit.named_colors()["black"], 
		linestyle = '--') 
	ax.scatter(radii, actual, marker = plots.mpltoolkit.markers()["x"], 
		color = plots.mpltoolkit.named_colors()["black"]) 
	bottom.plot(radii, len(radii) * [0], linestyle = '--', 
		c = plots.mpltoolkit.named_colors()["black"]) 
	bottom.scatter(radii, list(map(lambda x, y: x - y, actual, expected)), 
		color = plots.mpltoolkit.named_colors()["black"], 
		marker = plots.mpltoolkit.markers()["x"]) 


if __name__ == "__main__": 
	plt.clf() 
	ax, bottom = setup_axes() 
	out = vice.multioutput(sys.argv[1]) 
	# plot_expected(ax, out) 
	# plot_actual(ax, out) 
	plot_rates(ax, bottom, out) 
	plt.tight_layout() 
	plt.savefig(sys.argv[2]) 
	plt.clf() 

