r""" 
Produce a plot of :math:`\tau_\star` against time for each zone color-coded 
by radius. 

ARGV 
----
1) 		The name of the VICE output 
2) 		The name of the output image (without extension) 
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
XLIM = [-1, 14] 
YLIM = [0.2, 20] 

def setup_axis(): 
	r""" 
	Set up the matplotlib subplot to plot on 
	""" 
	fig = plt.figure(figsize = (7, 7)) 
	ax = fig.add_subplot(111, facecolor = "white") 
	ax.set_xlabel(r"Time [Gyr]") 
	ax.set_ylabel(r"$\tau_\star$ [Gyr]") 
	ax.set_xlim(XLIM) 
	ax.set_ylim(YLIM) 
	ax.set_yscale("log") 
	return ax 


def plot_tau_star(ax, out): 
	r""" 
	Plot :math:`\tau_\star` as a function of time in Gyr for each zone. 

	Parameters 
	----------
	ax : subplot 
		The matplotlib axes to plot on. 
	out : vice.multioutput 
		The multioutput object from the VICE simulation. 
	""" 
	cmap = plt.get_cmap(CMAP) 
	for i in range(int(15.5 / 0.25)): 
		tau_star = list(map(lambda x, y: 1.e-9 * x / y if y else float("nan"), 
			out.zones["zone%d" % (i)].history["mgas"], 
			out.zones["zone%d" % (i)].history["sfr"])) 
		ax.plot(out.zones["zone%d" % (i)].history["time"], 
			tau_star, c = cmap((0.25 * i + 0.125) / 15.5))
	sm = plt.cm.ScalarMappable(cmap = cmap, norm = plt.Normalize(vmin = 0, 
		vmax = 15.5)) 
	cbar = plt.colorbar(sm, cax = plots.mpltoolkit.append_axes(ax), pad = 0) 
	cbar.set_label(r"$R_\text{gal}$ [kpc]") 



if __name__ == "__main__": 
	plt.clf() 
	ax = setup_axis() 
	out = vice.output(sys.argv[1]) 
	plot_tau_star(ax, out) 
	plots.mpltoolkit.yticklabel_formatter(ax) 
	plt.tight_layout() 
	plt.savefig("%s.pdf" % (sys.argv[2])) 
	plt.savefig("%s.png" % (sys.argv[2])) 
	plt.clf() 

