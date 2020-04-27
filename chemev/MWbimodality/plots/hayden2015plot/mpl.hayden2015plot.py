r""" 
This script creates a heatmap in grayscale of the stellar number density in 
the [O/Fe]-[Fe/H] plane 

ARGV 
----
1) 		The name of the multizone output 
2) 		The name of the output image (without extension) 
""" 

import matplotlib.pyplot as plt 
import plots 
plots.mpltoolkit.load_mpl_presets() 
import numpy as np 
import vice 
import sys 

CMAP = "Greys" 
XLIM = [-1.2, 1.2] 
YLIM = [-0.3, 0.55] 


def setup_axes(): 
	r""" 
	Sets up the 3x5 matplotlib axes to plot on 
	""" 
	fig, axes = plt.subplots(ncols = 5, nrows = 3, figsize = (20, 12)) 
	for i in range(3): 
		for j in range(5): 
			if i != 2: plt.setp(axes[i][j].get_xticklabels(), visible = False) 
			if j != 0: plt.setp(axes[i][j].get_yticklabels(), visible = False) 
			axes[i][j].set_xlim(XLIM) 
			axes[i][j].set_ylim(YLIM) 
			if i == 0: axes[i][j].set_title(
				r"Final $R_\text{gal}$ = %g - %g kpc" % (
					[3, 5, 7, 9, 11][j], [5, 7, 9, 11, 13][j]), 
				fontsize = 25) 
			if j == 2: axes[i][j].text(-1, -0.2, 
				r"$\left|z\right|$ = %g - %g kpc" % (
					[1, 0.5, 0][i], [2, 1, 0.5][i]), 
				fontsize = 25) 
	axes[2][2].set_xlabel("[Fe/H]") 
	axes[1][0].set_ylabel("[O/Fe]") 
	return fig, axes 


def plot_density(ax, out, minrgal, maxrgal, minabsz, maxabsz): 
	r""" 
	Plot the 2D density of stars 

	ax : The matplotlib axis to plot the heatmap on 
	out : The VICE multioutput object 
	minrgal : Minimum galactocentric radius in kpc 
	maxrgal : Maximum galactocentric radius in kpc 
	minabsz : Minimum |z| in kpc 
	maxabsz : Maximum |z| in kpc 
	""" 
	stars = out.stars.filter("mass", ">", 0) 
	stars = stars.filter("zone_final", ">=", minrgal / 0.25) 
	stars = stars.filter("zone_final", "<=", (maxrgal - 0.25) / 0.25) 
	stars = stars.filter("abszfinal", ">=", minabsz) 
	stars = stars.filter("abszfinal", "<=", maxabsz) 
	dist, xedges, yedges = np.histogram2d(stars["[fe/h]"], stars["[o/fe]"], 
		range = [XLIM, YLIM], bins = 100, weights = stars["mass"], 
		density = True) 
	ax.imshow(np.rot90(dist), cmap = CMAP, aspect = "auto", 
		extent = [XLIM[0], XLIM[1], YLIM[0], YLIM[1]], interpolation = None, 
		vmin = 0, vmax = 15) 


if __name__ == "__main__": 
	plt.clf() 
	fig, axes = setup_axes() 
	out = vice.output(sys.argv[1]) 
	extra = np.genfromtxt("%s_extra_tracer_data.out" % (out.name)) 
	out.stars["abszfinal"] = [abs(row[-1]) for row in extra] 
	radii = [3, 5, 7, 9, 11, 13] 
	heights = [2, 1, 0.5, 0] 
	for i in range(3): 
		for j in range(5): 
			plot_density(axes[i][j], out, radii[j], radii[j + 1], 
				heights[i + 1], heights[i]) 
	# plot_density(axes[0][0], out, radii[0], radii[1], heights[1], heights[0]) 
	plt.tight_layout() 
	plt.subplots_adjust(hspace = 0, wspace = 0, left = 0.08) 
	plt.savefig("%s.pdf" % (sys.argv[2])) 
	plt.savefig("%s.png" % (sys.argv[2])) 
	plt.clf() 

