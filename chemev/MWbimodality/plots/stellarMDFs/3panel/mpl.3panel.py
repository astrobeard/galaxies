r""" 
Produces a plot of the stellar [Fe/H] distributions in 3 different z-bins 
and for different galactocentric radii. 

ARGV 
----
1) 		The name of the multizone output 
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

XLIM = [-1.4, 0.7] 
YLIM = [0, 0.19] 
BINS = np.linspace(-3, 1, 81).tolist() 


def setup_axes(): 
	r""" 
	Sets up the 3x1 array of matplotlib subplots to plot the distributions on 
	""" 
	fig = plt.figure(figsize = (7, 7)) 
	axes = 3 * [None] 
	for i in range(3): 
		axes[i] = fig.add_subplot(311 + i, facecolor = "white") 
		if i != 2: plt.setp(axes[i].get_xticklabels(), visible = False) 
		axes[i].set_xlim(XLIM) 
		axes[i].set_ylim(YLIM) 
		axes[i].yaxis.set_ticks([0, 0.05, 0.10, 0.15]) 
	axes[1].set_ylabel(r"PDF") 
	axes[-1].set_xlabel("[Fe/H]") 
	return axes 



def get_mdf(stars, minrgal, maxrgal, minabsz, maxabsz): 
	r""" 
	Calculates the stellar MDF in [Fe/H] for a given range in Rgal and |z| 

	Parameters 
	---------- 
	- stars : The star particles from the VICE output 
	- minrgal : The minimum galactocentric radius in kpc 
	- maxrgal : The maximum galactocentric radius in kpc 
	- minabsz : The minimum |z| in kpc 
	- maxabsz : The maximum |z| in kpc 
	""" 
	stars = stars.filter("zone_final", ">=", minrgal / 0.25) 
	stars = stars.filter("zone_final", "<=", (maxrgal - 0.25) / 0.25) 
	stars = stars.filter("abszfinal", ">=", minabsz) 
	stars = stars.filter("abszfinal", "<=", maxabsz) 
	stars = stars.filter("mass", ">", 0) 
	dist = (len(BINS) - 1) * [0.] 
	for i in range(len(dist)): 
		fltrd_stars = stars.filter("[fe/h]", ">=", BINS[i]) 
		fltrd_stars = fltrd_stars.filter("[fe/h]", "<=", BINS[i + 1]) 
		dist[i] = sum(fltrd_stars["mass"]) 
	norm = sum(dist) 
	return [i / norm for i in dist] 


def plot_mdf(ax, stars, minrgal, maxrgal, minabsz, maxabsz, color, 
	label = False): 
	r""" 
	Plot the MDF in [Fe/H] for a given range in Rgal and |z| 

	Parameters 
	----------
	- ax : The subplot to plot on 
	- stars : The star particles from the VICE output 
	- minrgal : The minimum galactocentric radius in kpc 
	- maxrgal : The maximum galactocentric radius in kpc 
	- minabsz : The minimum |z| in kpc 
	- maxabsz : The maximum |z| in kpc 
	- color : A string denoting the color to plot in 
	- label : Whether or not to label the lines for a legend 
	""" 
	centers = list(map(lambda x, y: (x + y) / 2, BINS[1:], BINS[:-1])) 
	mdf = get_mdf(stars, minrgal, maxrgal, minabsz, maxabsz) 
	kwargs = {} 
	if label: kwargs["label"] = r"%g $\leq R_\text{gal} \leq$ %g kpc" % (
		minrgal, maxrgal) 
	ax.plot(centers, mdf, c = plots.mpltoolkit.named_colors()[color], **kwargs)  


if __name__ == "__main__": 
	plt.clf() 
	axes = setup_axes() 
	out = vice.output(sys.argv[1]) 
	extra = np.genfromtxt("%s_extra_tracer_data.out" % (out.name)) 
	out.stars["abszfinal"] = [abs(row[-1]) for row in extra] 
	zbins = [2, 1, 0.5, 0] 
	rbins = [3, 5, 7, 9, 11, 13] 
	colors = ["red", "gold", "green", "blue", "darkviolet"] 
	for i in range(3): 
		axes[i].text(-0.6, 0.15, r"%g $\leq \left|z\right| \leq$ %g kpc" % (
			zbins[i + 1], zbins[i]), fontsize = 15)  
		for j in range(5): 
			plot_mdf(axes[i], out.stars, rbins[j], rbins[j + 1], 
				zbins[i + 1], zbins[i], colors[j], 
				label = i == 0) 
	leg = axes[0].legend(loc = plots.mpltoolkit.mpl_loc("upper left"), 
		ncol = 1, frameon = False, handlelength = 0, fontsize = 15) 
	for i in range(len(leg.get_texts())): 
		leg.get_texts()[i].set_color(colors[i]) 
		leg.legendHandles[i].set_visible(False) 
	plt.tight_layout() 
	plt.subplots_adjust(hspace = 0) 
	plt.savefig("%s.pdf" % (sys.argv[2])) 
	plt.savefig("%s.png" % (sys.argv[2])) 
	plt.clf() 

