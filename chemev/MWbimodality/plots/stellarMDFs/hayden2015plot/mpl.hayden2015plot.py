""" 
Plots the stellar [O/Fe] PDF in several [Fe/H] bins in a panel scheme similar 
to the Hayden et al. (2015) paper. 

ARGV 
----
1) 		The name of the multizone output 
2) 		The name of the output image (without the extension) 
""" 

import matplotlib.pyplot as plt 
import plots 
plots.mpltoolkit.load_mpl_presets() 
import numpy as np 
import math as m 
import vice 
import sys 
import os 

OFE_BINS = np.arange(-0.1, 0.51, 0.01).tolist() 
XLIM = [-0.12, 0.52] 
YLIM = [0.0, 0.24] 
FEH_BINS = [
	[-0.6, -0.4], 
	[-0.4, -0.2], 
	[-0.2, 0.0] 
] 
COLORS = [
	"dodgerblue", 
	"lime", 
	"crimson" 
]


def setup_axes(): 
	""" 
	Setup the 3x5 matplotlib axes. 
	""" 
	fig, axes = plt.subplots(ncols = 5, nrows = 3, figsize = (20, 12), 
		sharex = True, sharey = True)  
	for i in range(3): 
		for j in range(5): 
			if i != 2: plt.setp(axes[i][j].get_xticklabels(), visible = False) 
			if j != 0: plt.setp(axes[i][j].get_yticklabels(), visible = False) 
			axes[i][j].set_xlim(XLIM) 
			axes[i][j].set_ylim(YLIM) 
			axes[i][j].yaxis.set_ticks([0, 0.1, 0.2]) 
			axes[i][j].xaxis.set_ticks([0.0, 0.2, 0.4]) 
			if i == 0: axes[i][j].set_title(
				r"Final $R_\text{gal}$ = %g - %g kpc" % ( 
					[3, 5, 7, 9, 11][j], [5, 7, 9, 11, 13][j]), 
				fontsize = 25) 
			if j == 3: 
				axes[i][j].text(0.0, 0.2, 
					r"$\left|z\right|$ = %g - %g kpc" % ( 
						[1, 0.5, 0][i], [2, 1, 0.5][i]), 
					fontsize = 25) 
			elif j == 4: 
				axes[i][j].text(0.0, 0.2, 
					r"%g $\leq$ [Fe/H] $\leq$ %g" % (
						[-0.6, -0.4, -0.2][i], [-0.4, -0.2, 0.0][i]), 
					fontsize = 25, 
					color = plots.mpltoolkit.named_colors()[["crimson", "lime", 
						"dodgerblue"][i]])  
			else: pass 
	axes[2][2].set_xlabel("[O/Fe]") 
	axes[1][0].set_ylabel("PDF") 
	return axes 


def get_ofe_pdf(stars, min_rgal, max_rgal, minabsz, maxabsz, minFeH, maxFeH): 
	""" 
	Get the PDF within the binspace 

	stars : the dataframe of the stars from the VICE output 
	min_rgal : the lower bound galactocentric radius 
	max_rgal : the upper bound galactocentric radius 
	minabsz : the lower bound |z| 
	maxabsz : the upper bound |z| 
	min_FeH : The lower bound [Fe/H] to calculate the PDF for 
	max_FeH : The upper bound [Fe/H] to calculate the PDF for 
	"""	
	stars = stars.filter("zone_final", ">=", min_rgal / 0.25) 
	stars = stars.filter("zone_final", "<=", (max_rgal - 0.25) / 0.25) 
	stars = stars.filter("abszfinal", ">=", minabsz) 
	stars = stars.filter("abszfinal", "<=", maxabsz) 
	stars = stars.filter("[fe/h]", ">=", minFeH) 
	stars = stars.filter("[fe/h]", "<=", maxFeH) 
	if len(stars["mass"]) >= len(OFE_BINS): 
		dist = (len(OFE_BINS) - 1) * [0.0] 
		for i in range(len(dist)): 
			fltrd_stars = stars.filter("[o/fe]", ">=", OFE_BINS[i]) 
			fltrd_stars = fltrd_stars.filter("[o/fe]", "<=", OFE_BINS[i + 1]) 
			dist[i] += sum(fltrd_stars["mass"]) 
		norm = sum(dist) 
		return [i / norm for i in dist] 
	else: 
		return 


def plot_mdfs(ax, stars, min_rgal, max_rgal, minabsz, maxabsz, label = False): 
	""" 
	Plot all MDFs for a given rgal - |z| bin 

	ax : the subplot to plot on 
	stars : the dataframe of the stars from the VICE output 
	min_rgal : the minimum final galactocentric radius 
	max_rgal : the maximum final galactocentric radius 
	minabsz : the minimum final |z| 
	maxabsz : the maximum final |z| 
	""" 
	for i in range(len(FEH_BINS)): 
		dist = get_ofe_pdf(stars, min_rgal, max_rgal, minabsz, maxabsz, 
			FEH_BINS[i][0], FEH_BINS[i][1]) 
		if dist is not None: 
			kwargs = {
				"c": 		plots.mpltoolkit.named_colors()[COLORS[i]] 
			} 
			# if label: kwargs["label"] = r"%g $\leq$ [Fe/H] $\leq$ %g" % (
			if label: kwargs["label"] = r"[Fe/H] = %g - %g" % (
				FEH_BINS[i][0], FEH_BINS[i][1]) 
			ax.plot(
				list(map(lambda x, y: (x + y) / 2., OFE_BINS[1:], 
					OFE_BINS[:-1])), 
				dist, **kwargs) 
		else: 
			pass 


if __name__ == "__main__": 
	plt.clf() 
	axes = setup_axes() 
	out = vice.multioutput(sys.argv[1]) 
	extra_tracer_data = np.genfromtxt("%s_extra_tracer_data.out" % (out.name)) 
	out.stars["abszfinal"] = [abs(row[-1]) for row in 
		extra_tracer_data[:out.stars.size[0]]] 
	radii = [3, 5, 7, 9, 11, 13] 
	z = [2, 1, 0.5, 0] 
	for i in range(3): 
		for j in range(5): 
			# plot_mdfs(axes[i][j], out.stars, radii[j], radii[j + 1], 
			# 	z[i + 1], z[i], label = i == 2 and j == 4) 
			plot_mdfs(axes[i][j], out.stars, radii[j], radii[j + 1], 
				z[i + 1], z[i]) 
	leg = axes[0][4].legend(loc = plots.mpltoolkit.mpl_loc("upper right"), 
		ncol = 1, frameon = False, handlelength = 0, fontsize = 25) 
	for i in range(len(leg.get_texts())): 
		leg.get_texts()[i].set_color(COLORS[i]) 
		leg.legendHandles[i].set_visible(False) 
	plt.tight_layout() 
	plt.subplots_adjust(wspace = 0, hspace = 0) 
	plt.savefig("%s.pdf" % (sys.argv[2])) 
	plt.savefig("%s.png" % (sys.argv[2])) 
	plt.clf() 

