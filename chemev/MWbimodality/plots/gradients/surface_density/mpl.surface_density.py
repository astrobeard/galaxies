r""" 
Produce a plot of the surface density of gas and stars at the final timestep 
of a number of VICE simulations. 

ARGV 
----
1) 		The name of the output image (without extension) 
2-5) 	The name of the VICE outputs (up to 4) 
""" 

import matplotlib.pyplot as plt 
import plots 
plots.mpltoolkit.load_mpl_presets() 
import numpy as np 
import math as m 
import vice 
import sys 
import os 

RAD_BINS = np.linspace(0, 15.5, 63).tolist() 
COLORS = ["black", "crimson", "lime", "dodgerblue"] 
XLIM = [-1, 16] 
YLIM = [3e-4, 3e3] 


def setup_axis(): 
	fig = plt.figure(figsize = (7, 7)) 
	ax = fig.add_subplot(111, facecolor = "white") 
	ax.set_xlabel(r"$R_\text{gal}$ [kpc]") 
	ax.set_ylabel(r"Normalized $\Sigma$ [M$_\odot$ kpc$^{-2}$]") 
	ax.set_xlim(XLIM) 
	ax.set_ylim(YLIM) 
	ax.set_yscale("log") 
	return ax 


def plot_densities(ax, out, color): 
	stellar = (len(RAD_BINS) - 1) * [0.] 
	gaseous = (len(RAD_BINS) - 1) * [0.] 
	for i in range(len(stellar)): 
		stellar[i] = out.zones["zone%d" % (i)].history["mstar"][-1] 
		stellar[i] /= m.pi * (RAD_BINS[i + 1]**2 - RAD_BINS[i]**2) 
		gaseous[i] = out.zones["zone%d" % (i)].history["mgas"][-1] 
		gaseous[i] /= m.pi * (RAD_BINS[i + 1]**2 - RAD_BINS[i]**2) 
	norm_stellar = stellar[16] 
	norm_gaseous = gaseous[16] 
	for i in range(len(stellar)): 
		stellar[i] /= norm_stellar 
		gaseous[i] /= norm_gaseous 
	bin_centers = lambda x: [(x[i] + x[i + 1]) / 2 for i in range(len(x) - 1)] 
	ax.scatter(bin_centers(RAD_BINS), stellar, 
		marker = plots.mpltoolkit.markers()["star"], 
		c = plots.mpltoolkit.named_colors()[color], 
		s = 50) 
	ax.plot(bin_centers(RAD_BINS), gaseous, 
		linestyle = ':', 
		c = plots.mpltoolkit.named_colors()[color]) 


def legend(ax): 
	labels = (len(sys.argv) - 2) * [None] 
	for i in range(len(labels)): 
		if '/' in sys.argv[i + 2]: 
			labels[i] = sys.argv[i + 2].split('/')[-1].replace('_', r"\_") 
		else: 
			labels[i] = sys.argv[i + 2].replace('_', r"\_") 
	lines = len(labels) * [None] 
	for i in range(len(lines)): 
		lines[i] = ax.plot([1, 1], [2, 2], label = labels[i], 
			c = plots.mpltoolkit.named_colors()["white"])[0] 
	leg = ax.legend(loc = plots.mpltoolkit.mpl_loc("upper right"), ncol = 1, 
		frameon = False, bbox_to_anchor = (0.98, 0.98), handlelength = 0)  
	for i in range(len(lines)): 
		lines[i].remove() 
		leg.get_texts()[i].set_color(COLORS[i]) 


if __name__ == "__main__": 
	plt.clf() 
	ax = setup_axis() 
	for i in range(2, len(sys.argv)): 
		plot_densities(ax, vice.multioutput(sys.argv[i]), COLORS[i - 2]) 
	ax.plot(ax.get_xlim(), 
		[m.exp(-i / 3) * 1.e-1 for i in ax.get_xlim()], 
		c = plots.mpltoolkit.named_colors()["black"]) 
	legend(ax) 
	plots.mpltoolkit.yticklabel_formatter(ax) 
	plt.tight_layout() 
	plt.subplots_adjust(left = 0.2) 
	plt.savefig("%s.pdf" % (sys.argv[1])) 
	plt.savefig("%s.png" % (sys.argv[1])) 
	plt.clf() 

