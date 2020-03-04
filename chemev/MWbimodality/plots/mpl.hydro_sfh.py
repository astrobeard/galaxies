"""
This script produces a figure showing the star formation history that the 
hydro simulation had 

ARGV 
==== 
1) 		The name of the output image 
""" 

import matplotlib.pyplot as plt 
import matplotlib.colors as colors 
import matplotlib.cm as cm 
import plots 
plots.mpltoolkit.load_mpl_presets() 
import numpy as np 
import math as m 
import vice 
import sys 
import os 

XLIM = [-1, 15] 
CMAP = "plasma_r" 
RADIAL_BINS = np.linspace(0, 16, 17).tolist() 
TIME_BINS = np.linspace(0, 13.8, 36).tolist() 
DATA = np.genfromtxt("../data/UWhydro_particles.dat") 


def setup_axis(): 
	fig = plt.figure(figsize = (7, 7)) 
	ax = fig.add_subplot(111, facecolor = "white") 
	ax.set_xlabel("Time [Gyr]") 
	ax.set_ylabel(r"$\dot{M}_*$ [$M_\odot\ yr^{-1}$]") 
	ax.set_xlim(XLIM) 
	return ax 


def bin_centers(binspace): 
	return list(map(lambda x, y: (x + y) / 2., binspace[1:], binspace[:-1])) 


def plot_sfh_actual(ax, min_rgal, max_rgal, scalarMap): 
	sub = list(filter(lambda x: min_rgal <= x[2] <= max_rgal, DATA)) 
	counts, _ = np.histogram([i[1] for i in sub], bins = TIME_BINS) 
	counts = [i / 1000 for i in counts] 
	ax.plot(bin_centers(TIME_BINS), counts, 
		c = scalarMap.to_rgba((min_rgal + max_rgal) / 2)) 


def plot_sfh_measured(ax, min_rgal, max_rgal, scalarMap): 
	sub = list(filter(lambda x: min_rgal <= x[4] <= max_rgal, DATA)) 
	counts, _ = np.histogram([i[1] for i in sub], bins = TIME_BINS) 
	counts = [i / 1000 for i in counts] 
	ax.plot(bin_centers(TIME_BINS), counts, 
		c = scalarMap.to_rgba((min_rgal + max_rgal) / 2)) 


if __name__ == "__main__": 
	plt.clf() 
	ax = setup_axis() 
	cmap = plt.get_cmap(CMAP) 
	scalarMap = cm.ScalarMappable(norm = colors.Normalize(
		vmin = RADIAL_BINS[0], vmax = RADIAL_BINS[-1]), cmap = cmap) 
	for i in range(len(RADIAL_BINS) - 1): 
		plot_sfh_actual(ax, RADIAL_BINS[i], RADIAL_BINS[i + 1], scalarMap) 
		# plot_sfh_measured(ax, RADIAL_BINS[i], RADIAL_BINS[i + 1], scalarMap) 
	scalarMap._A = [] 
	cbar = plt.colorbar(scalarMap, ax = ax, pad = 0) 
	cbar.set_label(r"$R_\text{gal}$ [kpc]") 
	plt.tight_layout() 
	plt.savefig(sys.argv[1]) 
	plt.clf() 

