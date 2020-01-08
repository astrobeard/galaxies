""" 
Produces a plot of the APOGEE DR16 halo dwarfs and giants in the 
[alpha/M]-[M/H] plane. 

ARGV 
==== 
1)	The name of the output figure 
""" 

import matplotlib.pyplot as plt 
import plots 
plots.mpltoolkit.load_mpl_presets() 
import vice 
from vice.yields.presets import starburst19 
vice.yields.ccsne.settings["o"] = 0.01
import numpy as np 
import math as m 
import sys 
import os 

GIANTS_FILE = "../data/APOGEE/dr16_giants.dat" 
DWARFS_FILE = "../data/APOGEE/dr16_dwarfs.dat" 


def read_data(filename): 
	""" 
	Parameters 
	========== 
	filename :: str 
		The path to the file to read in 

	Returns 
	======= 
	The data in dictionary format 
	""" 
	raw = np.genfromtxt(filename, delimiter = ',')  
	data = {} 
	data["m_h"] 			= [row[6] for row in raw] 
	data["m_h_err"] 		= [row[7] for row in raw] 
	data["alpha_m"] 		= [row[8] for row in raw] 
	data["alpha_m_err"] 	= [row[9] for row in raw] 
	return data 


def setup_axis(): 
	""" 
	Creates and returns the matplotlib subplot object 
	""" 
	fig = plt.figure(figsize = (7, 7)) 
	ax = fig.add_subplot(111, facecolor = "white") 
	ax.set_xlabel(r"$\log_{10}(Z/Z_\odot)$") 
	ax.set_ylabel(r"[$\alpha$/M]") 
	ax.set_xlim([-2.8, 0]) 
	ax.set_ylim([-0.2, 0.8]) 
	return ax 


def plot_dataset(ax, data, color): 
	""" 
	Plots a dataset on the subplot 

	Parameters 
	========== 
	ax :: subplot 
		The matplotlib subplot object ot plot on 
	data :: dict 
		The data to plot from 
	marker :: str 
		The name of the matplotlib marker to plot the data in 
	color :: str 
		The name of the color to plot in 
	""" 
	ax.errorbar(data["m_h"], data["alpha_m"], 
		xerr = data["m_h_err"], yerr = data["alpha_m_err"], 
		c = plots.mpltoolkit.named_colors()[color], 
		fmt = 'o') 


def draw_legend(ax, colors, labels): 
	""" 
	Draws the legend on the axis 
	""" 
	lines = len(colors) * [None] 
	for i in range(len(lines)): 
		lines[i] = ax.plot([1, 2], [1, 2], 
			c = plots.mpltoolkit.named_colors()["white"], 
			label = labels[i])[0] 
	leg = ax.legend(loc = plots.mpltoolkit.mpl_loc("upper right"), ncol = 1, 
		frameon = False, bbox_to_anchor = (0.98, 0.98))  
	for i in leg.legendHandles: 
		i.set_visible(False) 
	for i in range(len(lines)): 
		lines[i].remove() 
		leg.get_texts()[i].set_color(colors[i])  


def run_vice_simulation(name = "onezonemodel", tau_star = 2, end = 1): 
	vice.singlezone(name = name, tau_star = tau_star, dt = 1e-3, eta = 1, 
		func = lambda t: 0).run(np.linspace(0, end, 1001), overwrite = True) 

def plot_vice_simulation(ax, name = "onezonemodel"): 
	hist = vice.history(name) 
	ax.plot(hist["[m/h]"], hist["[o/fe]"], 
		# list(map(lambda x, y: x - y, hist["[o/h]"], hist["[m/h]"])), 
		c = plots.mpltoolkit.named_colors()["black"], 
		zorder = 10) 

if __name__ == "__main__": 
	plt.clf() 
	giants = read_data(GIANTS_FILE) 
	dwarfs = read_data(DWARFS_FILE) 
	ax = setup_axis() 
	plot_dataset(ax, giants, "crimson") 
	plot_dataset(ax, dwarfs, "dodgerblue") 
	run_vice_simulation(name = "tau_star_2p5", tau_star = 2.5, end = 0.7) 
	run_vice_simulation(name = "tau_star_20", tau_star = 20, end = 2)  
	plot_vice_simulation(ax, "tau_star_2p5") 
	plot_vice_simulation(ax, "tau_star_20") 
	draw_legend(ax, ["crimson", "dodgerblue"], ["Giants", "Dwarfs"]) 
	plt.tight_layout()  
	plt.savefig(sys.argv[1]) 
	plt.clf() 


