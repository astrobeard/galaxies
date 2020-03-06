""" 
Plots the tracer particles from a given simulation on the [Y/X]-[X/H] 
axis for the 3-5, 5-7, 7-9, 9-11, and 11-13 kpc annuli within the disk 
color-coded by the age of each tracer particle. 

ARGV 
==== 
1) 	The name of the multizone output 
2)	The name of the output image 
3) 	The reference element X 
4) 	The secondary element Y 
""" 

import matplotlib.pyplot as plt 
import plots 
plots.mpltoolkit.load_mpl_presets() 
import numpy as np 
import math as m 
import vice 
import sys 
import os 

XLIM = [-1.2, 0.7] 
YLIM = [0.0, 0.5] 
CMAP = "plasma_r"

def setup_axes(): 
	fig = plt.figure(figsize = (35, 7)) 
	axes = 5 * [None] 
	for i in range(5): 
		axes[i] = fig.add_subplot(151 + i, facecolor = "white") 
		axes[i].set_title(r"%g kpc $\leq$ Final $R_\text{gal}$ $\leq$ %g kpc" % (
			[3, 5, 7, 9, 11][i], [5, 7, 9, 11, 13][i]), fontsize = 25) 
		axes[i].xaxis.set_ticks([-1.5, -1.0, -0.5, 0.0, 0.5]) 
		axes[i].set_xlim(XLIM) 
		axes[i].set_ylim(YLIM) 
		if i > 0: plt.setp(axes[i].get_yticklabels(), visible = False) 
	axes[0].set_ylabel("[%s/%s]" % (sys.argv[4], sys.argv[3])) 
	axes[2].set_xlabel("[%s/H]" % (sys.argv[3])) 
	return axes 

def plot_tracers(ax, tracers, zone_bounds): 
	cmap = plt.get_cmap(CMAP) 
	tracers = [list(i) for i in zip(
		tracers["formation_time"], 
		tracers["zone_final"], 
		tracers["mass"], 
		tracers["z(%s)" % (sys.argv[3])], 
		tracers["z(%s)" % (sys.argv[4])]
	)] 
	tracers = list(filter(lambda x: zone_bounds[0] <= x[1] <= zone_bounds[1], 
		tracers)) 
	tracers = list(filter(lambda x: x[3] > 0 and x[4] > 0, tracers)) 
	XH = len(tracers) * [0.] 
	YX = len(tracers) * [0.] 
	colors = len(tracers) * [None] 
	sizes = len(tracers) * [None] 
	for i in range(len(tracers)): 
		XH[i] = m.log10(tracers[i][3] / vice.solar_z[sys.argv[3]]) 
		YX[i] = m.log10(tracers[i][4] / vice.solar_z[sys.argv[4]]) - XH[i] 
		colors[i] = 13.8 - tracers[i][0] 
		sizes[i] = tracers[i][2] / 4e7 * 4 * (1 - 
			vice.cumulative_return_fraction(tracers[i][0])) 
	sc = ax.scatter(XH, YX, c = colors, s = sizes, cmap = cmap, vmin = 1, 
		vmax = 13.8) 
	return sc 


def plot_track(ax, szout): 
	ax.plot(szout.history["[%s/h]" % (sys.argv[3])], 
		szout.history["[%s/%s]" % (sys.argv[4], sys.argv[3])], 
		c = plots.mpltoolkit.named_colors()["lime"]) 


if __name__ == "__main__": 
	axes = setup_axes() 
	out = vice.multioutput(sys.argv[1]) 
	extra_tracer_data = np.genfromtxt("%s_extra_tracer_data.out" % (out.name)) 
	out.tracers["zfinal"] = [row[-1] for row in extra_tracer_data[:out.tracers.size[0]]]  
	fltrd_tracers = out.tracers.filter("zfinal", ">=", -3.) 
	fltrd_tracers = fltrd_tracers.filter("zfinal", "<=", 3.) 
	plot_tracers(axes[0], fltrd_tracers, [12, 19]) 
	plot_tracers(axes[1], fltrd_tracers, [20, 27]) 
	plot_tracers(axes[2], fltrd_tracers, [28, 35]) 
	plot_tracers(axes[3], fltrd_tracers, [36, 43]) 
	sc = plot_tracers(axes[4], fltrd_tracers, [44, 51]) 
	plot_track(axes[0], out.zones["zone12"]) 
	plot_track(axes[0], out.zones["zone19"]) 
	plot_track(axes[1], out.zones["zone20"]) 
	plot_track(axes[1], out.zones["zone27"])
	plot_track(axes[2], out.zones["zone28"]) 
	plot_track(axes[2], out.zones["zone35"])
	plot_track(axes[3], out.zones["zone36"]) 
	plot_track(axes[3], out.zones["zone43"])
	plot_track(axes[4], out.zones["zone44"]) 
	plot_track(axes[4], out.zones["zone51"])
	cbar = plt.colorbar(sc, ax = axes[4], pad = 0) 
	cbar.set_label("Age [Gyr]") 
	plt.tight_layout() 
	plt.subplots_adjust(wspace = 0) 
	plt.savefig(sys.argv[2], rasterized = True) 
	plt.clf() 

