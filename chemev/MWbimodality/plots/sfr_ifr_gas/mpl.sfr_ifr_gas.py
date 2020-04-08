r""" 
Produce a plot of the star formation rate, infall rate, and gas supply at all 
radii as a function of time for a given VICE output. 

ARGV 
----
1) 		The name of the VICE output 
2) 		The name of the output image (without extension) 
""" 

import matplotlib.pyplot as plt 
from matplotlib.ticker import FormatStrFormatter 
import plots 
plots.mpltoolkit.load_mpl_presets() 
import numpy as np 
import math as m 
import vice 
import sys 
import os 


CMAP = "plasma" 
XLIM = [-1, 14] 
RAD_BINS = np.linspace(0, 30, 121).tolist() 
YLOG = True 


def setup_axes(): 
	r""" 
	Set up the 1x3 matplotlib figure 
	""" 
	fig = plt.figure(figsize = (21, 7)) 
	axes = 3 * [None] 
	ylabels = [
		r"$\dot{M}_\text{in}$ [M$_\odot$ yr$^{-1}$]", 
		r"Normalized $\dot{M}_\star$ [M$_\odot$ yr$^{-1}$]", 
		r"Normalized $M_\text{gas}$ [M$_\odot$]"
	]
	for i in range(3): 
		axes[i] = fig.add_subplot(131 + i, facecolor = "white") 
		axes[i].set_xlabel(r"Time [Gyr]") 
		axes[i].set_ylabel(ylabels[i]) 
		axes[i].set_xlim(XLIM) 
		if YLOG: axes[i].set_yscale("log") 
	return [fig, axes] 


def plot_quantities(axes, out): 
	r""" 
	Plot the infall rate, star formation history, and gas mass against time in 
	the simulation. 

	Parameters 
	----------
	axes : list 
		The 3 axes to plot on. 
	out : vice.multioutput 
		The multioutput object from the VICE simulation. 
	""" 
	cmap = plt.get_cmap(CMAP) 
	# sigma_sfr_prefactor = 1 
	# sigma_gas_prefactor = 1e-9 
	for i in range(int(15.5 / 0.25)): 
		# sigma_sfr = out.zones["zone%d" % (i)].history.size[0] * [0.] 
		# sigma_gas = out.zones["zone%d" % (i)].history.size[0] * [0.] 
		# for j in range(out.zones["zone%d" % (i)].history.size[0]): 
		# 	area = m.pi * (RAD_BINS[i + 1]**2 - RAD_BINS[i]**2) 
		# 	sigma_sfr[j] = out.zones["zone%d" % (i)].history["sfr"][j] / area 
		# 	sigma_gas[j] = out.zones["zone%d" % (i)].history["mgas"][j] / area 
		# 	sigma_sfr[j] *= sigma_sfr_prefactor 
		# 	sigma_gas[j] *= sigma_gas_prefactor 

		kwargs = {
			"c": 		cmap((0.25 * i + 0.125) / 15.5) 
		} 
		axes[0].plot(out.zones["zone%d" % (i)].history["time"], 
			out.zones["zone%d" % (i)].history["ifr"], **kwargs) 
		axes[1].plot(out.zones["zone%d" % (i)].history["time"], 
			out.zones["zone%d" % (i)].history["sfr"], **kwargs) 
			# sigma_sfr, **kwargs) 
		axes[2].plot(out.zones["zone%d" % (i)].history["time"], 
			[i * 1.e-9 for i in out.zones["zone%d" % (i)].history["mgas"]], 
			**kwargs) 
			# sigma_gas, **kwargs) 


if __name__ == "__main__": 
	plt.clf() 
	fig, axes = setup_axes() 
	out = vice.output(sys.argv[1]) 
	plot_quantities(axes, out) 
	plt.tight_layout() 
	cbar_ax = fig.add_axes([0.92, 0.05, 0.02, 0.95]) 
	sm = plt.cm.ScalarMappable(cmap = plt.get_cmap(CMAP), 
		norm = plt.Normalize(vmin = 0, vmax = 15.5)) 
	cbar = plt.colorbar(sm, cax = cbar_ax) 
	cbar.set_label(r"$R_\text{gal}$ [kpc]") 
	# axes[1].yaxis.set_ticks(range(0, m.ceil(axes[1].get_ylim()[1]), 5)) 
	# axes[2].yaxis.set_ticks(range(0, m.ceil(axes[2].get_ylim()[1]), 2)) 
	axes[1].set_ylim([1e-3, 1e4]) 
	axes[2].set_ylim([1e-3, 1e4]) 
	plt.subplots_adjust(right = 0.92) 
	cbar_ax.set_position([
		axes[-1].get_position().x1, 
		axes[-1].get_position().y0, 
		0.02, 
		axes[-1].get_position().y1 - axes[-1].get_position().y0 
	]) 
	plt.savefig("%s.pdf" % (sys.argv[2])) 
	plt.savefig("%s.png" % (sys.argv[2])) 
	plt.clf() 

