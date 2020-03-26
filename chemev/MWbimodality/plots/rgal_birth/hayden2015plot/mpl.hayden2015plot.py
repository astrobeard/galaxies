""" 
Plots the tracer particles from a given simulation on the [Y/X]-[X/H] axis 
for the 3-5, 5-7, 7-9, 9-11, and 11-13 kps annuli within the disk and in 3 
rows for the |z| < 0.5, 0.5 < |z| < 1.0, and 1.0 < |z| < 2.0 vertical height 
bins. 

ARGV 
----
1) 		The name of the multizone output 
2) 		The name of the output figure (without the extension) 
""" 

import matplotlib.pyplot as plt 
import plots 
plots.mpltoolkit.load_mpl_presets() 
import numpy as np 
import math as m 
import vice 
import sys 
import os 

REF_ELEMENT = "Fe" 
SEC_ELEMENT = "O" 
XLIM = [-1.2, 0.7] 
YLIM = [-0.02, 0.52] 
CMAP = "plasma_r" 


def setup_axes(): 
	""" 
	Setup the 3x5 matplotlib axes. 
	""" 
	fig, axes = plt.subplots(ncols = 5, nrows = 3, figsize = (35, 21)) 
	for i in range(3): 
		for j in range(5): 
			if i != 2: plt.setp(axes[i][j].get_xticklabels(), visible = False) 
			if j != 0: plt.setp(axes[i][j].get_yticklabels(), visible = False) 
			axes[i][j].set_xlim(XLIM) 
			axes[i][j].set_ylim(YLIM) 
			if i == 0: axes[i][j].set_title(
				r"%g kpc $\leq$ Final $R_\text{gal}\ \leq$ %g kpc" % (
					[3, 5, 7, 9, 11][j], [5, 7, 9, 11, 13][j]), 
				fontsize = 25) 
			if j == 2: 
				axes[i][j].text(-0.5, 0.45, 
					r"%g kpc $\leq \left|z\right| \leq$ %g kpc" % (
						[1, 0.5, 0][i], [2, 1, 0.5][i]), 
					fontsize = 25) 
	axes[2][2].set_xlabel("[%s/H]" % (REF_ELEMENT)) 
	axes[1][0].set_ylabel("[%s/%s]" % (SEC_ELEMENT, REF_ELEMENT)) 
	return fig, axes 


def plot_stars(ax, stars, zone_bounds, zbounds): 
	""" 
	Plot the stars in a given radial range and z range on a given axis. 

	ax : subplot 
		The axes to plot on 
	stars : dataframe 
		The dataframe containing stars from the multioutput file (unfiltered) 
	zone_bounds : array-like 
		The inner and outer zones to take in the plot 
	zbounds : array-like 
		The lower and upper bound on |z| 
	""" 
	cmap = plt.get_cmap(CMAP) 
	stars = stars.filter("zone_final", ">=", zone_bounds[0]) 
	stars = stars.filter("zone_final", "<=", zone_bounds[1]) 
	stars = stars.filter("abszfinal", ">=", zbounds[0]) 
	stars = stars.filter("abszfinal", "<=", zbounds[1]) 
	colors = [i["zone_origin"] * 0.25 for i in stars] 
	sizes = [i["mass"] / 1e6 * 4 * (1 - 
		vice.cumulative_return_fraction(i["age"])) for i in stars] 
	return ax.scatter(
		stars["[%s/H]" % (REF_ELEMENT)], 
		stars["[%s/%s]" % (SEC_ELEMENT, REF_ELEMENT)], 
		c = colors, 
		s = sizes, 
		cmap = cmap, 
		vmin = 1, 
		vmax = 15 
	)  


if __name__ == "__main__": 
	plt.clf() 
	fig, axes = setup_axes() 
	out = vice.multioutput(sys.argv[1]) 
	extra_tracer_data = np.genfromtxt("%s_extra_tracer_data.out" % (out.name)) 
	out.stars["abszfinal"] = [abs(row[-1]) for row in 
		extra_tracer_data[:out.stars.size[0]]] 
	zone_bounds = [[12, 19], [20, 27], [28, 35], [36, 43], [44, 51]] 
	z_bounds = [[1, 2], [0.5, 1], [0, 0.5]] 
	for i in range(3): 
		for j in range(5): 
			sc = plot_stars(axes[i][j], out.stars, zone_bounds[j], z_bounds[i]) 
	cbar_ax = fig.add_axes([0.95, 0.05, 0.02, 0.95]) 
	fig.colorbar(sc, cax = cbar_ax) 
	cbar_ax.set_ylabel(r"$R_\text{gal}$ of birth [kpc]") 
	plt.tight_layout() 
	plt.subplots_adjust(wspace = 0, hspace = 0, right = 0.95) 
	cbar_ax.set_position([
		axes[-1][-1].get_position().x1, 
		axes[-1][-1].get_position().y0, 
		0.02, 
		axes[0][-1].get_position().y1 - axes[-1][-1].get_position().y0
	])
	plt.savefig("%s.pdf" % (sys.argv[2])) 
	plt.savefig("%s.png" % (sys.argv[2])) 
	plt.clf() 
