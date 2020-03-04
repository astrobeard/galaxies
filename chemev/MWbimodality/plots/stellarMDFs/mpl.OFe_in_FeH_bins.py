""" 
This script produces a plot of the stellar [O/Fe] distribution in several 
[Fe/H] bins for a given disk model a radial bin in each panel 

ARGV 
==== 
1) 		The name of the multizone output 
2) 		The name of the output image 
""" 

import matplotlib.pyplot as plt 
from matplotlib.ticker import FormatStrFormatter as fsf 
import plots 
plots.mpltoolkit.load_mpl_presets() 
import numpy as np 
import math as m 
import vice 
import sys 
import os 

OFE_BINS = np.arange(0., 0.51, 0.01).tolist() 
XLIM = [-0.02, 0.52] 
YLIM = [0.0, 0.25] 
FEH_BINS = [ 
	# [-0.8, -0.7], 
	[-0.6, -0.5], 
	[-0.4, -0.3], 
	[-0.2, -0.1] 
] 
COLORS = [
	# "black", 
	"dodgerblue", 
	"lime", 
	"crimson" 
]


def setup_axes(): 
	""" 
	Setup the subplots 
	""" 
	fig = plt.figure(figsize = (35, 7)) 
	axes = 5 * [None] 
	for i in range(5): 
		axes[i] = fig.add_subplot(151 + i, facecolor = "white") 
		axes[i].set_title(r"%g kpc $\leq$ Final $R_\text{gal}$ $\leq$ %g kpc" % (
			[3, 5, 7, 9, 11][i], [5, 7, 9, 11, 13][i]), fontsize = 25) 
		# axes[i].set_yscale("log") 
		axes[i].set_xlim(XLIM) 
		axes[i].set_ylim(YLIM) 
		if i > 0: plt.setp(axes[i].get_yticklabels(), visible = False) 
	axes[0].yaxis.set_major_formatter(fsf("%g"))  
	axes[0].set_ylabel("Stellar PDF")  
	axes[2].set_xlabel("[O/Fe]") 
	return axes 


def get_ofe_pdf(stars, min_rgal, max_rgal, min_FeH, max_FeH): 
	""" 
	stars :: The dataframe of the stars from the VICE output 
	min_rgal :: The lower bound galactocentric radius 
	max_rgal :: The upper bound galactocentric radius 
	min_FeH :: The lower bound [Fe/H] to calculate the PDF for 
	max_FeH :: The upper bound [Fe/H] to calculate the PDF for 
	""" 
	stars = stars.filter("zone_final", ">=", min_rgal / 0.25) 
	stars = stars.filter("zone_final", "<=", (max_rgal - 0.25) / 0.25) 
	stars = stars.filter("[fe/h]", ">=", min_FeH) 
	stars = stars.filter("[fe/h]", "<=", max_FeH) 
	if len(stars["mass"]) >= len(OFE_BINS): 
		dist = (len(OFE_BINS) - 1) * [0.] 
		for i in range(len(dist)): 
			fltrd_stars = stars.filter("[o/fe]", ">=", OFE_BINS[i]) 
			fltrd_stars = fltrd_stars.filter("[o/fe]", "<=", OFE_BINS[i + 1]) 
			dist[i] += sum(fltrd_stars["mass"]) 
		norm = sum(dist) 
		return [i / norm for i in dist] 
	else: 
		return None 


def plot_pdfs(ax, stars, min_rgal, max_rgal, label = False): 
	for i in range(len(FEH_BINS)): 
		dist = get_ofe_pdf(stars, min_rgal, max_rgal, FEH_BINS[i][0], 
			FEH_BINS[i][1]) 
		if dist is not None: 
			kwargs = {
				"c": 	plots.mpltoolkit.named_colors()[COLORS[i]] 
			} 
			if label: kwargs["label"] = r"%g $\leq$ [Fe/H] $\leq$ %g" % (
				FEH_BINS[i][0], FEH_BINS[i][1]) 
			ax.plot(
				list(map(lambda x, y: (x + y) / 2., OFE_BINS[1:], OFE_BINS[:-1])), 
				dist, **kwargs) 
		else: 
			pass 


if __name__ == "__main__": 
	plt.clf() 
	axes = setup_axes() 
	out = vice.multioutput(sys.argv[1]) 
	extra_tracer_data = np.genfromtxt("%s_extra_tracer_data.out" % (out.name)) 
	for i in extra_tracer_data: 
		if i[-1] == 100: i[-1] = 0 
	out.tracers["zfinal"] = [row[-1] for row in extra_tracer_data[:out.tracers.size[0]]] 
	stars = out.tracers.filter("zfinal", ">=", -3.) 
	stars = stars.filter("zfinal", "<=", 3.) 
	stars = stars.filter("mass", ">", 0) 
	print("Number of stars: %d" % (len(stars["mass"]))) 
	plot_pdfs(axes[0], stars, 3, 5) 
	plot_pdfs(axes[1], stars, 5, 7) 
	plot_pdfs(axes[2], stars, 7, 9) 
	plot_pdfs(axes[3], stars, 9, 11) 
	plot_pdfs(axes[4], stars, 11, 13, label = True) 
	leg = axes[4].legend(loc = plots.mpltoolkit.mpl_loc("upper right"), ncol = 1, 
		frameon = False, handlelength = 0, fontsize = 25) 
	for i in range(len(leg.get_texts())): 
		leg.get_texts()[i].set_color(COLORS[i]) 
	plt.tight_layout() 
	plt.subplots_adjust(wspace = 0) 
	plt.savefig(sys.argv[2]) 
	plt.clf() 

