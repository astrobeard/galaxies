""" 
This script produces a figure showing PDFs of the birth radius of stars for 
different age ranges and in bins of final radius. 

ARGV 
==== 
1) 		The name of the output figure 
""" 

from importlib import import_module 
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
formation_bins = np.linspace(0, 16, 51).tolist() 
centers = list(map(lambda x, y: (x + y) / 2, formation_bins[1:], 
	formation_bins[:-1]))  

def setup_axes(): 
	fig, axes = plt.subplots(nrows = 3, ncols = 5, 
		sharex = True, sharey = True, figsize = (35, 21)) 
	axes[1][0].set_ylabel(r"PDF of birth $R_\text{gal}$", fontsize = 40) 
	axes[2][2].set_xlabel(r"$R_\text{gal}$ [kpc]", fontsize = 40) 
	for i in range(3): 
		for j in range(5): 
			axes[i][j].text(0.2, 1.0, 
				r"Final $R_\text{gal}$ %d - %d kpc" % (5 * i + j, 
					5 * i + j + 1), 
				fontsize = 40) 
	axes[0][0].set_ylim([0.0, 1.2]) 
	return axes 


def get_formation_radii(age_range, radius_range): 
	of_interest = len(UWhydro["rfinal"]) * [0] 
	for i in range(len(UWhydro["rfinal"])): 
		if radius_range[0] <= UWhydro["rfinal"][i] <= radius_range[1]: 
			if age_range[0] <= 13.8 - UWhydro["tform"][i] <= age_range[1]: 
				of_interest[i] = 1 
			else: pass 
		else: pass 
	formation = sum(of_interest) * [None] 
	n = 0 
	for i in range(len(UWhydro["rfinal"])): 
		if of_interest[i]: 
			formation[n] = UWhydro["rform"][i] 
			n += 1 
		else: continue 
	return formation 


def plot_formation_radii_pdfs(ax, radius_range, age_binspace, colors): 
	for i in range(len(age_binspace) - 1): 
		rforms = get_formation_radii([age_binspace[i], age_binspace[i + 1]], 
			radius_range) 
		pdf, bins = np.histogram(rforms, bins = formation_bins, 
			density = True) 
		ax.plot(centers, pdf, c = plots.mpltoolkit.named_colors()[colors[i]]) 
	inside = 0. 
	outside = 0. 
	insitu = 0. 
	rforms = get_formation_radii([0, 14], radius_range) 
	inside = len(list(filter(lambda x: x <= radius_range[0], rforms))) 
	outside = len(list(filter(lambda x: x >= radius_range[1], rforms))) 
	insitu = len(list(filter(lambda x: radius_range[0] <= x <= radius_range[1], 
		rforms)))  
	s = inside + outside + insitu 
	inside /= s 
	outside /= s 
	insitu /= s 
	ax.text(0.2, 0.9, r"from inside: %.2f" % (inside), fontsize = 30) 
	ax.text(0.2, 0.8, r"from outside: %.2f" % (outside), fontsize = 30) 
	ax.text(0.2, 0.7, r"from in-situ: %.2f" % (insitu), fontsize = 30) 
	ax.plot(2 * [radius_range[0]], ax.get_ylim(), 
		c = plots.mpltoolkit.named_colors()["black"], 
		linestyle = ':') 
	ax.plot(2 * [radius_range[1]], ax.get_ylim(), 
		c = plots.mpltoolkit.named_colors()["black"], 
		linestyle = ':') 



if __name__ == "__main__": 
	plt.clf() 
	axes = setup_axes() 
	age_binspace = [0, 2, 4, 6, 8, 10, 12, 14] 
	colors = ["black", "grey", "darkviolet", "dodgerblue", "lime", "gold", 
		"crimson"] 
	for i in range(2): 
		for j in range(3): 
			axes[i][j].text(8, 0.6, "Age = %d - %d Gyr" % (
				age_binspace[3 * i + j], age_binspace[3 * i + j + 1]), 
			color = plots.mpltoolkit.named_colors()[colors[3 * i + j]], 
			fontsize = 30) 
	axes[2][0].text(8, 0.6, "Age = 12 - 14 Gyr", 
		color = plots.mpltoolkit.named_colors()["crimson"], 
		fontsize = 30) 
	for i in range(3): 
		for j in range(5): 
			plot_formation_radii_pdfs(axes[i][j], 
				[5 * i + j, 5 * i + j + 1], 
				age_binspace, colors)  
	plt.tight_layout() 
	plt.subplots_adjust(wspace = 0, hspace = 0)  
	plt.savefig(sys.argv[1]) 
	plt.clf() 

