""" 
Plots the tracer particles from a given simulation on the [Y/X]-[X/H] 
axis 

ARGV 
==== 
1) 	The name of the VICE output, without the '.vice' extension 
2) 	The name of the output image 
3) 	The reference element X 
4) 	The secondary element Y 
5) 	The column number of Z(X) in the output tracer particle file 
6) 	The column number of Z(Y) in the output tracer particle file 
""" 


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
	ax.set_xlabel("[%s/H]" % (sys.argv[3])) 
	ax.set_ylabel("[%s/%s]" % (sys.argv[4], sys.argv[3])) 
	return ax 


def tracer_data(): 
	return np.genfromtxt("%s.vice/tracers.out" % (sys.argv[1])) 


def plot_tracers(ax, tracers): 
	cmap = plt.get_cmap("viridis") 
	tracers = list(filter(lambda x: 40 <= x[2] <= 55, tracers)) 
	tracers = list(filter(lambda x: x[5] > 0, tracers)) 
	tracers = list(filter(lambda x: x[6] > 0, tracers)) 
	XH = len(tracers) * [0.] 
	YX = len(tracers) * [0.] 
	colors = len(tracers) * [0.] 
	sizes = len(tracers) * [0.] 
	for i in range(len(tracers)): 
		XH[i] = m.log10(tracers[i][int(sys.argv[5])] / 
			vice.solar_z[sys.argv[3]]) 
		YX[i] = m.log10(tracers[i][int(sys.argv[6])] / 
			vice.solar_z[sys.argv[4]]) - XH[i] 
		colors[i] = tracers[i][1] * 0.25 
		sizes[i] = 20 * tracers[i][3] / 4e6 
		sys.stdout.write("Progress: %.2f%%\r" % (100. * (i + 1) / len(tracers))) 
		sys.stdout.flush() 
	sys.stdout.write("\n") 
	sc = ax.scatter(XH, YX, c = colors, s = sizes, cmap = cmap) 
	cbar = plt.colorbar(sc, ax = ax, pad = 0) 
	cbar.set_label(r"$R_\text{gal}$ of birth [kpc]") 
	ax.set_title(r"10 kpc $\leq$ Final $R_\text{gal}$ $\leq$ 14 kpc", 
		fontsize = 25) 


if __name__ == "__main__": 
	plt.clf() 
	ax = setup_axis() 
	plot_tracers(ax, tracer_data()) 
	plt.tight_layout() 
	ax.set_xlim([-1.7, 0.2]) 
	ax.set_ylim([-0.24, 0.24]) 
	plt.savefig(sys.argv[2]) 
	plt.clf() 



