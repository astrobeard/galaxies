""" 
This script produces a figure showing PDFs of the final height above the plane 
and the vertical velocity distributions. 

ARGV 
==== 
1) 			The name of the output figure 
""" 

import matplotlib.pyplot as plt 
import plots 
plots.mpltoolkit.load_mpl_presets() 
import numpy as np 
import math as m 
import sys 
import os 
current = os.getcwd() 
os.chdir("../../") 
sys.path.append(os.getcwd()) 
os.chdir(current) 
from data import UWhydro 

zfinal_bins = np.linspace(-10, 10, 1001).tolist() 
vz_bins = np.linspace(-300, 300, 1001).tolist() 

def setup_axes(): 
	fig = plt.figure(figsize = (14, 7)) 
	axes = 2 * [None] 
	for i in range(len(axes)): 
		axes[i] = fig.add_subplot(121 + i, facecolor = "white") 
		axes[i].set_ylabel("PDF") 
		axes[i].set_yscale("log") 
	axes[0].set_xlabel(r"$z_\text{final}$ [kpc]") 
	axes[1].set_xlabel(r"$v_\text{z}$ [km s$^{-1}$]") 
	axes[0].set_xlim([-12, 12]) 
	axes[1].set_xlim([-350, 350]) 
	return axes 


def bin_centers(arr): 
	return list(map(lambda x, y: (x + y) / 2., arr[1:], arr[:-1])) 


def plot_zfinal_dist(ax): 
	pdf, _ = np.histogram(UWhydro["zfinal"], bins = zfinal_bins, density = True) 
	ax.plot(bin_centers(zfinal_bins), pdf, 
		c = plots.mpltoolkit.named_colors()["black"]) 
	ax.set_ylim(ax.get_ylim()) # fix the y-axis limits 
	ax.plot(2 * [-3], ax.get_ylim(), linestyle = ':', 
		c = plots.mpltoolkit.named_colors()["crimson"]) 
	ax.plot(2 * [3], ax.get_ylim(), linestyle = ':', 
		c = plots.mpltoolkit.named_colors()["crimson"]) 


def plot_vz_dist(ax): 
	pdf, _ = np.histogram(UWhydro["v_z"], bins = vz_bins, density = True) 
	ax.plot(bin_centers(vz_bins), pdf, 
		c = plots.mpltoolkit.named_colors()["black"]) 
	ax.set_ylim(ax.get_ylim()) # fix the y-axis limits 
	ax.plot(2 * [-50], ax.get_ylim(), linestyle = ':', 
		c = plots.mpltoolkit.named_colors()["crimson"]) 
	ax.plot(2 * [50], ax.get_ylim(), linestyle = ':', 
		c = plots.mpltoolkit.named_colors()["crimson"]) 


if __name__ == "__main__": 
	plt.clf() 
	axes = setup_axes() 
	plot_zfinal_dist(axes[0]) 
	plot_vz_dist(axes[1]) 
	plt.tight_layout() 
	plt.savefig(sys.argv[1]) 
	plt.clf() 


