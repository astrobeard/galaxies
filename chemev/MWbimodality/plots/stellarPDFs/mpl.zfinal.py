""" 
This script produces PDFs of zfinal for stars that migrated inward, outward, 
and in-situ for several ages in a given radial bin. 

ARGV 
==== 
1) 		The name of the output image 
2)		zfinal or vz, for which quantity to plot the PDF for 
3) 		The minimum galactocentric radius in kpc 
4)		The maximum galactocentric radius in kpc 
""" 

import matplotlib.pyplot as plt 
import plots 
plots.mpltoolkit.load_mpl_presets() 
import numpy as np 
import sys 
import os 
if os.getcwd() != os.path.dirname(os.path.abspath(__file__)): 
	os.chdir(os.path.dirname(os.path.abspath(__file__))) 
sys.path.append("../../") 
from data import UWhydro 

AGE_BINS = [[2, 3], [6, 7], [10, 11]] 
BINS = np.linspace(-250, 250, 101).tolist() 
XLIM = [-260, 260]


def setup_axes(): 
	""" 
	Setup the subplots 
	""" 
	fig = plt.figure(figsize = (len(AGE_BINS) * 7, 7)) 
	axes = len(AGE_BINS) * [None] 
	for i in range(len(axes)): 
		kwargs = {
			"facecolor": 		"white" 
		} 
		if i > 0: kwargs["sharey"] = axes[0] 
		axes[i] = fig.add_subplot(101 + 10 * len(AGE_BINS) + i, **kwargs) 
		if i > 0: plt.setp(axes[i].get_yticklabels(), visible = False) 
		axes[i].set_title(r"%g Gyr $\leq$ age $\leq$ %g Gyr" % (
			AGE_BINS[i][0], AGE_BINS[i][1]), fontsize = 25) 
		axes[i].set_xlim(XLIM) 
	axes[0].set_ylabel("Stellar PDF") 
	if sys.argv[2].lower() in ["zfinal", "v_z"]: 
		axes[len(AGE_BINS) // 2].set_xlabel({
			"zfinal":		r"$z_\text{final}$ [kpc]", 
			"v_z": 			r"$v_z$ [km s$^{-1}$]" 	
		}[sys.argv[2].lower()]) 
	else: 
		raise ValueError("ARGV[2] must be either 'zfinal' or 'v_z'. Got: %s" % (
			sys.argv[2])) 
	return axes 


def get_zfinal_pdf(stars): 
	""" 
	stars :: The VICE dataframe holding the stars in a given radial and age bin 
	""" 
	return np.histogram(stars[sys.argv[2].lower()], bins = BINS, 
		density = True)[0] 


def bin_centers(binspace): 
	""" 
	binspace :: an array of bin edges 
	""" 
	return list(map(lambda x, y: (x + y) / 2, binspace[1:], binspace[:-1])) 


def plot_zfinal_pdf(ax, stars, label = False): 
	insitu = stars.filter("rfinal", ">=", float(sys.argv[3])) 
	insitu = insitu.filter("rfinal", "<=", float(sys.argv[4])) 
	movedin = stars.filter("rfinal", "<=", float(sys.argv[3])) 
	movedout = stars.filter("rfinal", ">=", float(sys.argv[4])) 
	kwargs = {
		"c":		"black" 
	} 
	if label: kwargs["label"] = "in-situ" 
	ax.plot(bin_centers(BINS), get_zfinal_pdf(insitu), **kwargs) 
	kwargs["c"] = "crimson" 
	if label: kwargs["label"] = "moved in" 
	ax.plot(bin_centers(BINS), get_zfinal_pdf(movedin), **kwargs) 
	kwargs["c"] = "dodgerblue" 
	if label: kwargs["label"] = "moved out" 
	ax.plot(bin_centers(BINS), get_zfinal_pdf(movedout), **kwargs) 


if __name__ == "__main__": 
	plt.clf() 
	axes = setup_axes() 
	fltrd = UWhydro.filter("rform", ">=", float(sys.argv[3])) 
	fltrd = fltrd.filter("rform", "<=", float(sys.argv[4])) 
	stars = len(AGE_BINS) * [None] 
	for i in range(len(stars)): 
		stars[i] = fltrd.filter("tform", "<=", 12.8 - AGE_BINS[i][0]) 
		stars[i] = stars[i].filter("tform", ">=", 12.8 - AGE_BINS[i][1]) 
		plot_zfinal_pdf(axes[i], stars[i], label = i == 0) 
	leg = axes[0].legend(loc = plots.mpltoolkit.mpl_loc("upper left"), ncol = 1, 
		bbox_to_anchor = (0.01, 0.99), frameon = False, handlelength = 0) 
	for i in range(3): 
		leg.get_texts()[i].set_color(["black", "crimson", "dodgerblue"][i]) 
	plt.tight_layout() 
	plt.subplots_adjust(wspace = 0) 
	plt.savefig(sys.argv[1]) 
	plt.clf() 
