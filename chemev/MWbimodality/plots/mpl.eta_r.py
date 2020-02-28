""" 
This plot produces a figure of eta(r) as the original linear trend and 
subsequently the mode of [Mg/H] trend 

ARGV 
==== 
1) 		The name of the output image 
""" 

import matplotlib.pyplot as plt 
import plots 
plots.mpltoolkit.load_mpl_presets() 
import numpy as np 
import vice 
import sys 

# The yield of magnesium assuming our O yield and [O/Mg] = 0 
Y_MG_CC = 0.015 * vice.solar_z["mg"] / vice.solar_z["o"] 


def setup_axis(): 
	""" 
	Sets up the matplotlib subplot to plot on 
	""" 
	fig = plt.figure(figsize = (7, 7)) 
	ax = fig.add_subplot(111, facecolor = "white") 
	ax.set_xlabel(r"$R_\text{gal}$ [kpc]") 
	ax.set_ylabel(r"$\eta(r) = \dot{M}_\text{out} / \dot{M}_*$") 
	ax.set_xlim([-2, 32]) 
	ax.set_ylim([-2, 32]) 
	return ax 


def linear(rgal): 
	""" 
	The linear trend of eta with galactocentric radius 

	rgal :: The galactocentric radius in kpc 
	""" 
	return 0.8 * rgal 


def mode_mgh(rgal): 
	""" 
	The mode of [Mg/H] motivated eta with galactocentric radius 

	rgal :: The galactocentric radius in kpc 
	""" 
	return Y_MG_CC / vice.solar_z["mg"] * 10**(0.06 * (rgal - 4) - 0.3) - 0.6 


def plot_trends(ax): 
	""" 
	Plot the eta(r) trends 

	ax :: The subplot to plot on 
	""" 
	radii = np.linspace(0, 30, 1000) 
	ax.plot(radii, list(map(linear, radii)), 
		c = plots.mpltoolkit.named_colors()["crimson"]) 
	ax.plot(radii, list(map(mode_mgh, radii)), 
		c = plots.mpltoolkit.named_colors()["dodgerblue"]) 


if __name__ == "__main__": 
	plt.clf() 
	ax = setup_axis() 
	plot_trends(ax) 
	plt.tight_layout() 
	plt.savefig(sys.argv[1]) 
	plt.clf() 

