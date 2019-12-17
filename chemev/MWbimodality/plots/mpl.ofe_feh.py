""" 
Plots the tracer particles from a given simulation on the [O/Fe]-[Fe/H] 
axis 

ARGV 
==== 
1) 	The name of the VICE output, without the '.vice' extension 
2) 	The name of the output image 
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
	ax.set_xlabel("[Fe/H]") 
	ax.set_ylabel("[O/Fe]") 
	return ax 


def tracer_data(): 
	return np.genfromtxt("%s.vice/tracers.out" % (sys.argv[1])) 


def plot_tracers(ax, tracers): 
	cmap = plt.get_cmap("viridis") 
	tracers = list(filter(lambda x: x[2] in [7, 8], tracers)) 
	tracers = list(filter(lambda x: x[5] > 0, tracers)) 
	tracers = list(filter(lambda x: x[6] > 0, tracers)) 
	for i in range(len(tracers)): 
		# if tracers[i][2] in [7, 8]: 
		# 	try: 
		FeH = m.log10(tracers[i][5] / vice.solar_z["fe"]) 
		OH = m.log10(tracers[i][6] / vice.solar_z["o"]) 
		ax.scatter(FeH, OH - FeH, c = cmap(tracers[i][1] / 30), 
			s = 20 * tracers[3] / 4e6) 
		# 	except ValueError: 
		# 		continue 
		# else: 
		# 	pass 
		sys.stdout.write("Progress: %.2f%%\r" % (100. * (i + 1) / len(tracers))) 
		sys.stdout.flush() 
	sys.stdout.write("\n") 


if __name__ == "__main__": 
	plt.clf() 
	ax = setup_axis() 
	plot_tracers(ax, tracer_data()) 
	plt.tight_layout() 
	ax.set_xlim([-1.7, 0.2]) 
	ax.set_ylim([0., 0.5]) 
	plt.savefig(sys.argv[2]) 
	plt.clf() 



