""" 
This script creates a plot showing a heatmap in polar coordinates of a given 
quantity from a multizone simulation. Galactocentric radius is shown as the 
radial coordinate with time as the angular coordinate. 

ARGV 
==== 
1) 		The name of the output image 
2)		The name of the multizone output 
""" 

from vice.yields.presets import my_yields 
import matplotlib 
matplotlib.use("Agg") 
import matplotlib.pyplot as plt 
import matplotlib as mpl 
import plots 
plots.mpltoolkit.load_mpl_presets() 
import numpy as np 
import vice 
import sys 
import os 

CMAP = "bwr" 
# KEY = "[o/fe]" 
LABEL = r"\% difference $\dot{N}_\text{Ia}$"  
VMIN = -50 
VMAX = 50 


def get_heatmap(out): 
	radii = [0.25 * i for i in range(62)] 
	times = out.zones["zone0"].history["time"][:] 
	qty = len(radii) * [None] 
	for i in range(len(qty)): 
		# qty[i] = out.zones["zone%d" % (i)].history[KEY][:] 
		actual = get_proxies(out.zones["zone%d" % (i)]) 
		sz = vice.mirror(out.zones["zone%d" % (i)]) 
		sz.func = lambda t: out.zones["zone%d" % (i)].history["mgas"][0] 
		sz.name = "comparison" 
		comp = sz.run(np.linspace(0, 12.8, 641), overwrite = True, 
			capture = True) 
		expected = get_proxies(comp) 
		qty[i] = list(map(lambda x, y: 100 * (x - y) / y if y > 0 else 0, 
			actual, expected)) 
	return [radii, times, qty] 


# From the Ia rate calculation functions 
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
			proxies[i] -= (
				vice.yields.ccsne.settings["fe"] * zone.history["sfr"][i] * 1e9 
			) 
			proxies[i] += zone.history["z(fe)"][i] * zone.history["sfr"][i] * (
				1 + mir.eta - 0.4) * 1e9 
			# proxies[i] /= zone.history["mass(fe)"][i] 
			if proxies[i] < 0: proxies[i] = 0
		else: pass 
	return proxies 


def setup_axis(): 
	""" 
	Sets up the polar axis 
	""" 
	fig = plt.figure(figsize = (7, 7)) 
	ax = fig.add_subplot(111, facecolor = "white") 
	ax.set_xlabel("Time [Gyr]") 
	ax.set_ylabel(r"$R_\text{gal}$ [kpc]") 
	ax.yaxis.set_ticks([0, 3, 6, 9, 12, 15]) 
	return ax 


def draw_heatmap(ax, radii, times, qty): 
	im = ax.imshow(qty, cmap = CMAP, aspect = "auto", origin = "lower", 
		extent = [times[0], times[-1], radii[0], radii[-1]], 
		vmin = VMIN, vmax = VMAX, interpolation = None) 
	cbar = plt.gcf().colorbar(im, ax = ax, pad = 0.08) 
	cbar.set_label(LABEL)


if __name__ == "__main__": 
	plt.clf() 
	ax = setup_axis() 
	out = vice.output(sys.argv[2]) 
	radii, times, qty = get_heatmap(out) 
	draw_heatmap(ax, radii, times, qty)  
	plt.tight_layout() 
	plt.savefig(sys.argv[1]) 
	plt.clf() 


