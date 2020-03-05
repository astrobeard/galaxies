""" 
This script produces a plot of the Ia rate in several annular zones as a 
function of time in comparison to what is expected given the star formation 
history 

ARGV 
==== 
1) 		The name of the multizone output 
2) 		The name of the output image 
""" 

from vice.yields.presets import my_yields 
import matplotlib.pyplot as plt 
import plots 
plots.mpltoolkit.load_mpl_presets() 
import numpy as np 
import math as m 
import vice 
import sys 
import os 


def setup_axis(): 
	fig = plt.figure(figsize = (7, 7)) 
	ax = fig.add_subplot(111, facecolor = "white") 
	ax.set_ylabel(r"$m_\text{Fe}^\text{Ia}\dot{N}_\text{Ia}/M_\text{Fe}$ [Gyr$^{-1}$]") 
	ax.set_xlabel("Time [Gyr]") 
	ax.set_xlim([-1, 15]) 
	return ax 


def get_proxies(zone): 
	mir = vice.mirror(zone) 
	proxies = (len(zone.history["time"]) - 1) * [0.] 
	for i in range(len(proxies)): 
		if zone.history["time"][i] > mir.delay: 
			proxies[i] = (
				zone.history["mass(fe)"][i + 1] - zone.history["mass(fe)"][i] 
			) / (
				zone.history["time"][i + 1] - zone.history["time"][i] 
			) 
			if i == 8: print(proxies[i]) 
			proxies[i] -= (
				vice.yields.ccsne.settings["fe"] * zone.history["sfr"][i] * 1e9 
			) 
			if i == 8: print(proxies[i]) 
			proxies[i] += zone.history["z(fe)"][i] * zone.history["sfr"][i] * (
				1 + mir.eta - 0.4) * 1e9 
			if i == 8: print(proxies[i]) 
			proxies[i] /= zone.history["mass(fe)"][i] 
			if i == 8: print(proxies[i]) 
			if proxies[i] < 0: proxies[i] = 0
		else: pass 
	return proxies 


def plot_actual(ax, zone, color, norm, label): 	
	# proxies = [i / norm for i in get_proxies(zone)] 
	proxies = get_proxies(zone) 
	kwargs = {
		"c": 		plots.mpltoolkit.named_colors()[color] 
	}
	if label is not None: kwargs["label"] = label 
	print(proxies[:15]) 
	ax.plot(zone.history["time"][:-1], proxies, **kwargs) 


def plot_comparison(ax, zone, color): 
	sz = vice.mirror(zone) 
	sz.func = lambda t: zone.history["mgas"][0] 
	sz.name = "comparison" 
	comp = sz.run(np.linspace(0, 14, 1401), overwrite = True, capture = True) 
	# proxies = [i / comp.history["mass(fe)"][-1] for i in get_proxies(comp)] 
	proxies = get_proxies(comp) 
	print(proxies[:15]) 
	ax.plot(comp.history["time"][:-1], proxies, 
		c = plots.mpltoolkit.named_colors()[color], linestyle = '--') 
	return comp.history["mass(fe)"][-1] 



if __name__ == "__main__": 
	plt.clf() 
	ax = setup_axis() 
	out = vice.output(sys.argv[1]) 
	radii = [5, 10, 15]
	colors = ["dodgerblue", "lime", "crimson"] 
	for i in range(len(radii)): 
		norm = plot_comparison(ax, out.zones["zone%d" % (int(radii[i] / 0.25))], 
			colors[i]) 
		plot_actual(ax, out.zones["zone%d" % (int(radii[i] / 0.25))], 
			colors[i], norm, r"$R_\text{gal}$ = %g kpc" % (radii[i])) 
	leg = ax.legend(loc = plots.mpltoolkit.mpl_loc("lower right"), ncol = 1, 
		frameon = False, bbox_to_anchor = (0.99, 0.01), handlelength = 0) 
	for i in range(len(leg.legendHandles)): 
		leg.get_texts()[i].set_color(colors[i]) 
	plt.tight_layout() 
	plt.savefig(sys.argv[2]) 
	plt.clf() 




