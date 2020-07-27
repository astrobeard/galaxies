""" 
This script produces a plot of the Ia rate in several annular zones as a 
function of time in comparison to what is expected given the star formation 
history 

ARGV 
==== 
1) 		The name of the multizone output 
2) 		The name of the output image (without extension) 
""" 

from vice.yields.presets import JW20  
import matplotlib.pyplot as plt 
import plots 
plots.mpltoolkit.load_mpl_presets() 
import numpy as np 
import math as m 
import vice 
import sys 
import os 
sys.path.append("../../simulations/") 
import gas_disks 
from conference import TIME_SWITCH, tau_in, tau_star, eta 


def setup_axis(): 
	fig = plt.figure(figsize = (7, 7)) 
	ax = fig.add_subplot(111, facecolor = "white") 
	ax.set_ylabel(
		r"$m_\text{Fe}^\text{Ia}\dot{N}_\text{Ia}/M_\text{Fe}$ [Gyr$^{-1}$]") 
	ax.set_xlabel("Time [Gyr]") 
	ax.set_xlim([-1, 15]) 
	ax.set_ylim([-0.1, 1.3])
	return ax 


def get_proxies(zone): 
	mir = vice.singlezone.from_output(zone) 
	proxies = (len(zone.history["time"]) - 1) * [0.] 
	for i in range(len(proxies)): 
		if zone.history["time"][i] > mir.delay: 
			proxies[i] = (
				zone.history["mass(fe)"][i + 1] - zone.history["mass(fe)"][i] 
			) / (
				zone.history["time"][i + 1] - zone.history["time"][i] 
			) 
			proxies[i] -= (
				vice.yields.ccsne.settings["fe"] * zone.history["sfr"][i] * 1e9 
			) 
			proxies[i] += zone.history["z(fe)"][i] * zone.history["sfr"][i] * (
				1 + mir.eta - 0.4) * 1e9 
			proxies[i] /= zone.history["mass(fe)"][i] 
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
	ax.plot(zone.history["time"][:-1], proxies, **kwargs) 


def plot_comparison(ax, zone, color): 
	sz = vice.singlezone.from_output(zone) 
	# sz.func = lambda t: zone.history["mgas"][0] 
	# sz.schmidt = False 
	# sz.func = vice.toolkit.repair_function(zone, "sfr") 
	# sz.tau_star = vice.toolkit.repair_function(zone, "tau_star") 
	sz.name = "comparison" 
	comp = sz.run(np.linspace(0, 14, 1401), overwrite = True, capture = True) 
	# proxies = [i / comp.history["mass(fe)"][-1] for i in get_proxies(comp)] 
	proxies = get_proxies(comp) 
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
		norm = plot_comparison(ax, 
			out.zones["zone%d" % (int(radii[i] / 0.25))], colors[i]) 
		# norm = 1 
		plot_actual(ax, out.zones["zone%d" % (int(radii[i] / 0.25))], 
			colors[i], norm, r"$R_\text{gal}$ = %g kpc" % (radii[i])) 
	leg = ax.legend(loc = plots.mpltoolkit.mpl_loc("lower right"), ncol = 1, 
		frameon = False, bbox_to_anchor = (0.99, 0.01), handlelength = 0) 
	for i in range(len(leg.legendHandles)): 
		leg.get_texts()[i].set_color(colors[i]) 
		leg.legendHandles[i].set_visible(False) 
	plt.tight_layout() 
	plt.savefig("%s.pdf" % (sys.argv[2])) 
	plt.savefig("%s.png" % (sys.argv[2])) 
	plt.clf() 




