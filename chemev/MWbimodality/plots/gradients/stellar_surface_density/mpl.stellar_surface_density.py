""" 
This script produces a plot of the stellar surface density as a function of 
radius for a given VICE multizone output. 

ARGV 
---- 
1) 		The name of the multizone output 
2) 		The name of the output image (no extension) 
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

RAD_BINS = np.linspace(0, 30, 121).tolist() 
NORM = 1e11 


def setup_axis(): 
	""" 
	Sets up a subplot to plot the stellar surface density on 
	""" 
	fig = plt.figure(figsize = (7, 7)) 
	ax = fig.add_subplot(111, facecolor = "white") 
	ax.set_xlabel(r"$R_\text{gal}$ [kpc]") 
	ax.set_ylabel(r"$\propto\Sigma_\star$ [M$_\odot$ pc$^{-2}$]") 
	ax.set_xlim([-1, 21]) 
	ax.set_ylim([0.0003, 300]) 
	ax.set_yscale("log") 
	ax.yaxis.set_major_formatter(fsf("%g")) 
	return ax 


def surface_density(output): 
	""" 
	Compute the stellar surface densities 

	output : vice.multioutput 
		The multioutput object with stellar data 
	""" 
	densities = (len(RAD_BINS) - 1) * [0.] 
	stars = output.stars.filter("mass", ">", 0) 
	for i in range(len(stars["mass"])): 
		densities[int(stars["zone_final"][i])] += stars["mass"][i] * (1 - 
			vice.cumulative_return_fraction(stars["age"][i])) 
	for i in range(len(densities)): 
		densities[i] /= m.pi * (RAD_BINS[i + 1]**2 - RAD_BINS[i]**2) 
		densities[i] /= NORM 
	return densities 


def draw(ax, output): 
	""" 
	Plot the stellar surface densities as a function of radius 

	ax : subplot 
		The subplot to plot on 
	output : vice.multioutput 
		The multioutput object with stellar data 
	""" 
	bin_centers = lambda x: list(map(lambda x, y: (x + y) / 2., x[1:], x[:-1])) 
	ax.scatter(bin_centers(RAD_BINS), surface_density(output), 
		marker = plots.mpltoolkit.markers()["star"], 
		c = plots.mpltoolkit.named_colors()["black"]) 
	xvals = np.linspace(0, 21, 1001).tolist() 
	Rs = 2 
	profile = lambda x: 13 * m.exp(-x / Rs) 
	ax.plot(xvals, [profile(i) for i in xvals], linestyle = ':', 
		c = plots.mpltoolkit.named_colors()["black"]) 
	ax.text(11, 20, r"$R_{\text{s},\star} \approx$ 2 kpc", fontsize = 25) 
	ax.text(11, 50, r"$R_\text{s,g}$ = 3 kpc", fontsize = 25)  


if __name__ == "__main__": 
	plt.clf() 
	ax = setup_axis() 
	out = vice.multioutput(sys.argv[1]) 
	draw(ax, out) 
	plt.tight_layout() 
	plt.savefig("%s.pdf" % (sys.argv[2])) 
	plt.savefig("%s.png" % (sys.argv[2])) 
	plt.clf() 

