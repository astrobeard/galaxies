""" 
Produces a plot comparing model tracks and predicited [alpha/M] distributions 
to the APOGEE+Gaia observations. 

ARGV 
==== 
1)	The name of the output figure 
2) 	1 to rerun the VICE simulations, 0 if not necessary 
""" 

import matplotlib.pyplot as plt 
import plots 
plots.mpltoolkit.load_mpl_presets() 
import vice 
from vice.yields.presets import starburst19 
vice.yields.ccsne.settings["o"] = 0.0115
from data import apogee_gaia 
import numpy as np 
import math as m 
import sys 
import os 

def setup_axes(): 
	fig = plt.figure(figsize = (14, 7)) 
	axes = 2 * [None] 
	for i in range(2): 
		axes[i] = fig.add_subplot(121 + i, facecolor = "white") 
	axes[0].set_xlabel(r"$\log_{10}(Z/Z_\odot)$") 
	axes[0].set_ylabel(r"[$\alpha$/M]") 
	axes[1].set_xlabel(r"[$\alpha$/M]") 
	axes[1].set_ylabel(r"PDF") 
	# axes[1].set_yscale("log") 
	axes[0].set_xlim([-2.8, 0]) 
	axes[0].set_ylim([-0.2, 0.8]) 
	axes[1].set_xlim([-0.25, 0.65]) 
	return axes 

def scatter_stars(ax, data): 
	ax.errorbar(data["m_h"], data["alpha_m"], 
		xerr = data["m_h_err"], yerr = data["alpha_m_err"], 
		c = plots.mpltoolkit.named_colors()["black"], 
		marker = '.', ms = 1, linestyle = '')  

def dist_stars(ax, data): 
	normed, bins = np.histogram(data["alpha_m"], 
		bins = np.linspace(-0.2, 0.8, 21), density = True) 
	ax.step(bins[1:], normed, c = plots.mpltoolkit.named_colors()["black"])

# def run_vice_simulation(name = "onezonemodel", tau_star = 2, end = 1): 
# 	vice.singlezone(name = name, tau_star = tau_star, dt = 1e-4, eta = 1,  
# 		bins = np.linspace(-3, 1, 201), 
# 		func = lambda t: 0).run(np.linspace(0, end, 1001), overwrite = True) 

def run_vice_simulation(name = "onezonemodel", tau_in = 1, end = 1, 
	tau_star = 2): 
	kwargs = {} 
	kwargs["name"] 		= name 
	kwargs["tau_star"] 	= tau_star
	kwargs["schmidt"] 	= True 
	kwargs["dt"] 		= 1e-4
	kwargs["eta"] 		= 1
	kwargs["bins"] 		= np.linspace(-3, 1, 401) 
	kwargs["Mg0"] 		= 0
	kwargs["func"]		= lambda t: 6.0 * m.exp(-t / tau_in) 
	kwargs["MgSchmidt"] = 6.0e9
	sz = vice.singlezone(**kwargs) 
	print(sz) 
	sz.run(np.linspace(0, end, 1001), overwrite = True) 


def plot_vice_track(ax, name = "onezonemodel", linestyle = '-', 
	color = "black"): 
	out = vice.output(name) 
	ax.plot(out.history["[fe/h]"], 
		out.history["[o/fe]"], 
		c = plots.mpltoolkit.named_colors()[color], 
		zorder = 10) 

def plot_vice_distributions(ax, names, proportions, colors): 
	outputs = [vice.output(i) for i in names]
	for i in range(len(names)): 
		centers = list(map(lambda x, y: (x + y) / 2, 
			outputs[i].mdf["bin_edge_left"], outputs[i].mdf["bin_edge_right"])) 
		ax.plot(centers, outputs[i].mdf["dn/d[o/fe]"], 
			c = plots.mpltoolkit.named_colors()[colors[i]], 
			zorder = 10) 
	total = len(outputs[0].mdf["bin_edge_left"]) * [0.] 
	for i in range(len(outputs[0].mdf["bin_edge_left"])): 
		for j in range(len(names)): 
			total[i] += outputs[j].mdf["dn/d[o/fe]"][i] * proportions[j] 
	convolved = convolve_total(centers, total, 0.03) 
	ax.plot(list(map(lambda x, y: (x + y) / 2, 
		outputs[0].mdf["bin_edge_left"], outputs[0].mdf["bin_edge_right"])), 
		total, c = plots.mpltoolkit.named_colors()["black"], 
		zorder = 12, linestyle = ':') 
	ax.plot(list(map(lambda x, y: (x + y) / 2, 
		outputs[0].mdf["bin_edge_left"], outputs[0].mdf["bin_edge_right"])), 
		convolved, c = plots.mpltoolkit.named_colors()["black"], 
		zorder = 12, linestyle = '--') 

def convolve_total(centers, dist, disp): 
	convolved = len(centers) * [0.] 
	for i in range(len(centers)): 
		s = 0
		for j in range(len(centers)): 
			s += dist[j] * m.exp(-((centers[j] - centers[i]) / disp)**2 / 2) 
		convolved[i] = s 

	s = 0 
	for i in range(len(convolved)): 
		s += convolved[i] * (centers[1] - centers[0]) 
	for i in range(len(convolved)): 
		convolved[i] /= s 
	return convolved 

if __name__ == "__main__": 
	plt.clf() 
	data = apogee_gaia.whole() 
	axes = setup_axes() 
	scatter_stars(axes[0], data) 
	dist_stars(axes[1], data) 
	names = ["tau_in_0p5", "tau_in_0p2"] 
	timescales = [0.5, 0.2] 
	if int(sys.argv[2]): 
		# run_vice_simulation(name = "tau_star_1", tau_star = 1, end = 0.4) 
		# run_vice_simulation(name = "tau_star_8", tau_star = 8, end = 2) 
		run_vice_simulation(name = names[0], tau_in = timescales[0], end = 1.5, 
			tau_star = 2.5) 
		run_vice_simulation(name = names[1], tau_in = timescales[1], end = 3, 
			tau_star = 15) 
	else: pass 
	plot_vice_track(axes[0], name = names[0], color = "crimson") 
	plot_vice_track(axes[0], name = names[1], color = "dodgerblue")   
	plot_vice_distributions(axes[1], [names[0], names[1]], 
		[0.5, 0.5], ["crimson", "dodgerblue"]) 
	plt.tight_layout() 
	plt.savefig(sys.argv[1]) 
	plt.clf() 


