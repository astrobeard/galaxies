
import matplotlib.pyplot as plt 
from matplotlib.ticker import FormatStrFormatter as fsf 
import plots 
plots.mpltoolkit.load_mpl_presets() 
from astropy.io import fits 
import numpy as np 
import vice 
import sys 

OUTPUTSDIR = "/Users/astrobeard/Work/Research/VICErepos/VICE/migration/outputs" 
STATIC = "%s/high-resolution/2Gyr/diffusion/static" % (OUTPUTSDIR)
INSIDEOUT = "%s/high-resolution/2Gyr/diffusion/insideout" % (OUTPUTSDIR) 
LATEBURST = "%s/high-resolution/2Gyr/diffusion/lateburst" % (OUTPUTSDIR) 
OUTERBURST = "%s/high-resolution/2Gyr/diffusion/outerburst" % (OUTPUTSDIR)
# INSIDEOUT = "../../simulations/paper_noburst" 
# LATEBURST = "../../simulations/paper_withburst" 
STEM = "age_metallicity" 
CMAP = "winter" 
ZONE_WIDTH = 0.1 
ZONE_MIN = int(7 / ZONE_WIDTH) 
ZONE_MAX = int((9 - ZONE_WIDTH) / ZONE_WIDTH) 
LOGAGE = True 


def setup_axes(): 
	fig = plt.figure(figsize = (20, 10)) 
	axes = 2 * [None] 
	for i in range(len(axes)): 
		axes[i] = 4 * [None] 
		for j in range(len(axes[i])): 
			axes[i][j] = fig.add_subplot(241 + 4 * i + j, facecolor = "white") 
			if i != 1: plt.setp(axes[i][j].get_xticklabels(), visible = False) 
			if j != 0: plt.setp(axes[i][j].get_yticklabels(), visible = False) 
			if LOGAGE: 
				axes[i][j].set_xscale("log") 
				axes[i][j].set_xlim([0.2, 18]) 
				axes[i][j].xaxis.set_major_formatter(fsf("%g")) 
			else: 
				axes[i][j].set_xlim([-1, 14]) 
			if i: 
				axes[i][j].set_ylim([-1.1, 0.7]) 
				axes[i][j].yaxis.set_ticks([-1. + 0.2 * i for i in range(9)]) 
			else: 
				axes[i][j].set_ylim([-0.9, 0.7]) 
				axes[i][j].yaxis.set_ticks([-0.8 + 0.2 * i for i in range(8)]) 
				axes[i][j].text(0.5, -0.7, 
					["Static", "Inside-Out", "Late-Burst", "Outer-Burst"][j], 
					fontsize = 25) 
		axes[i][0].set_ylabel(["[O/H]", "[Fe/H]"][i]) 
	dummy = fig.add_subplot(111, facecolor = "white", zorder = -1) 
	posd = dummy.get_position() 
	posd.x0 = axes[1][0].get_position().x0 
	posd.x1 = axes[1][2].get_position().x1 
	posd.y0 = axes[1][0].get_position().y0 
	posd.y1 = axes[0][0].get_position().y1 
	dummy.set_position(posd) 
	dummy.set_xlabel("Age [Gyr]", labelpad = 30) 
	plt.setp(dummy.get_xticklabels(), visible = False) 
	plt.setp(dummy.get_yticklabels(), visible = False) 
	axes.append(dummy) 

	return axes 


def plot_amr(ax, element, output): 
	cmap = plt.get_cmap(CMAP) 
	stars = output.stars.filter("zone_final", ">=", ZONE_MIN) 
	stars = stars.filter("zone_final", "<=", ZONE_MAX) 
	stars = stars.filter("zfinal", ">=", -0.5) 
	stars = stars.filter("zfinal", "<=", 0.5) 
	stars = stars.filter("mass", ">=", 1.) 
	colors = [ZONE_WIDTH * (i + 0.5) for i in stars["zone_origin"]] 
	return ax.scatter(stars["age"], stars["[%s/H]" % (element)], c = colors, 
		s = 0.1, cmap = cmap, vmin = 0, vmax = 15) 


def zheights(name): 
	raw = np.genfromtxt("%s_analogdata.out" % (name)) 
	return [row[-1] for row in raw] 


def weighted_median(values, weights, stop = 0.5): 
	indeces = np.argsort(values) 
	values = [values[i] for i in indeces] 
	weights = [weights[i] for i in indeces] 
	weights = [i / sum(weights) for i in weights] 
	s = 0 
	for i in range(len(weights)): 
		s += weights[i] 
		if s > stop: 
			idx = i - 1 
			break 
	return values[idx] 


def median_ages(ax, element, output, label = False): 
	stars = output.stars.filter("zone_final", ">=", ZONE_MIN) 
	stars = stars.filter("zone_final", "<=", ZONE_MAX) 
	stars = stars.filter("zfinal", ">=", -0.5) 
	stars = stars.filter("zfinal", "<=", 0.5) 
	stars = stars.filter("mass", ">=", 1.) 
	bins = np.arange(-1., 1.05, 0.05) 
	ages = (len(bins) - 1) * [0.] 
	lowers = (len(bins) - 1) * [0.] 
	uppers = (len(bins) - 1) * [0.] 
	for i in range(len(ages)): 
		stars_ = stars.filter("[%s/h]" % (element), ">=", bins[i]) 
		stars_ = stars_.filter("[%s/h]" % (element), "<=", bins[i + 1]) 
		if len(stars_["age"]) > 20: 
			masses = list(map(lambda x, y: x * (1 - 
				vice.cumulative_return_fraction(y)), 
			stars_["mass"], stars_["age"])) 
			ages[i] = weighted_median(stars_["age"], masses) 
			# ages[i] = np.mean(stars_["age"]) 
			lowers[i] = weighted_median(stars_["age"], masses, stop = 0.16) 
			uppers[i] = weighted_median(stars_["age"], masses, stop = 0.84) 
		else: 
			ages[i] = float("nan") 
			lowers[i] = float("nan") 
			uppers[i] = float("nan") 
	# ax.scatter(ages, list(map(lambda x, y: (x + y) / 2., bins[1:], bins[:-1])), 
	# 	marker = plots.mpltoolkit.markers()["star"], 
	# 	c = plots.mpltoolkit.named_colors()["black"], s = 100) 
	xerr = [ 
		[ages[i] - lowers[i] for i in range(len(ages))], 
		[uppers[i] - ages[i] for i in range(len(ages))] 
	] 
	kwargs = { 
		"xerr": 		xerr, 
		"yerr": 		(bins[1] - bins[0]) / 2, 
		"c": 			plots.mpltoolkit.named_colors()["black"], 
		"marker": 		plots.mpltoolkit.markers()["square"], 
		"linestyle": 	"None" 
	} 
	if label: kwargs["label"] = "Simulation" 
	ax.errorbar(ages, list(map(lambda x, y: (x + y) / 2, bins[1:], bins[:-1])), 
		**kwargs) 


def feuillet2018_data(ax, element, label = False): 
	raw = np.genfromtxt({
			"fe": 		"age_mh.dat", 
			"o": 		"age_oh.dat" 
		}[element.lower()]) 
	onh = len(raw) * [0.] 
	onh_disp = len(raw) * [0.] 
	age = len(raw) * [0.] 
	age_disp = [len(raw) * [0.], len(raw) * [0.]] 
	for i in range(len(raw)): 
		onh[i] = (raw[i][0] + raw[i][1]) / 2 
		onh_disp[i] = (raw[i][1] - raw[i][0]) / 2 
		age[i] = 10**(raw[i][2] - 9) 
		age_disp[0][i] = age[i] - 10**(raw[i][2] - raw[i][3] - 9) 
		age_disp[1][i] = 10**(raw[i][2] + raw[i][3] - 9) - age[i] 
	# ax.scatter(age, onh, marker = plots.mpltoolkit.markers()["triangle_up"], 
	# 	c = plots.mpltoolkit.named_colors()["crimson"], s = 100) 
	kwargs = {
		"xerr": 		age_disp, 
		"yerr": 		onh_disp, 
		"c": 			plots.mpltoolkit.named_colors()["crimson"], 
		"marker": 		plots.mpltoolkit.markers()["triangle_up"], 
		"linestyle": 	"None" 
	} 
	if label: kwargs["label"] = "Feuillet et al. (2018)" 
	ax.errorbar(age, onh, **kwargs) 


def feuillet2019_data(ax, element, label = False): 
	raw = fits.open({
		"fe": 	"Feuillet2019_MH_ages/ELEM_GAUSS_AGE_07_09_00_05_M_H.fits", 
		"o": 	"APOGEE_DR14_OH_ages/ELEM_GAUSS_AGE_07_09_00_05_O_H.fits" 
	}[element.lower()]) 
	onh = len(raw[1].data) * [0.] 
	onh_disp = len(raw[1].data) * [0.] 
	age = len(raw[1].data) * [0.] 
	age_disp = [len(raw[1].data) * [0.], len(raw[1].data) * [0.]] 
	for i in range(len(raw[1].data)): 
		if raw[1].data['nstars'][i] > 15: 
			onh[i] = (raw[1].data['bin_ab'][i] + 
				raw[1].data['bin_ab_max'][i]) / 2 
			onh_disp[i] = (raw[1].data['bin_ab_max'][i] - 
				raw[1].data['bin_ab'][i]) / 2 
			age[i] = 10**(raw[1].data['mean_age'][i] - 9) 
			age_disp[0][i] = age[i] - 10**(raw[1].data['mean_age'][i] - 
				raw[1].data['age_disp'][i] - 9) 
			age_disp[1][i] = 10**(raw[1].data['mean_age'][i] + 
				raw[1].data['age_disp'][i] - 9) - age[i] 
		else: 
			onh[i] = float("nan") 
			onh_disp[i] = float("nan") 
			age[i] = float("nan") 
			age_disp[0][i] = float("nan") 
			age_disp[1][i] = float("nan") 
	kwargs = {
		"xerr": 		age_disp, 
		"yerr": 		onh_disp, 
		"c": 			plots.mpltoolkit.named_colors()["crimson"], 
		"marker": 		plots.mpltoolkit.markers()["triangle_up"], 
		"linestyle": 	"None" 
	} 
	if label: kwargs["label"] = "Feuillet et al. (2019)" 
	ax.errorbar(age, onh, **kwargs) 


if __name__ == "__main__": 
	axes = setup_axes() 
	static = vice.output(STATIC) 
	insideout = vice.output(INSIDEOUT) 
	lateburst = vice.output(LATEBURST) 
	outerburst = vice.output(OUTERBURST) 
	static.stars["zfinal"] = zheights(STATIC)[:static.stars.size[0]] 
	insideout.stars["zfinal"] = zheights(INSIDEOUT)[:insideout.stars.size[0]] 
	lateburst.stars["zfinal"] = zheights(LATEBURST)[:lateburst.stars.size[0]] 
	outerburst.stars["zfinal"] = zheights(OUTERBURST)[:outerburst.stars.size[0]] 
	for i in range(2): 
		for j in range(4): 
			sc = plot_amr(axes[i][j], ["O", "Fe"][i], 
				[static, insideout, lateburst, outerburst][j]) 
			median_ages(axes[i][j], ["O", "Fe"][i], 
				[static, insideout, lateburst, outerburst][j], 
				label = i == 1 and j == 0) 
			feuillet2019_data(axes[i][j], ["O", "Fe"][i], 
				label = i == 1 and j == 0) 

	axes[1][0].legend(loc = plots.mpltoolkit.mpl_loc("lower left"), ncol = 1, 
		frameon = False, bbox_to_anchor = (0.01, 0.01), fontsize = 18) 
	cbar_ax = plt.gcf().add_axes([0.92, 0.05, 0.02, 0.95]) 
	cbar = plt.colorbar(sc, cax = cbar_ax, pad = 0.0, orientation = "vertical") 
	cbar.set_label(r"$R_\text{gal}$ of birth [kpc]", labelpad = 10) 
	cbar.set_ticks(range(2, 16, 2)) 
	plt.tight_layout() 
	plt.subplots_adjust(wspace = 0, hspace = 0, bottom = 0.1, right = 0.9) 
	cbar_ax.set_position([
		axes[1][-1].get_position().x1, 
		axes[1][-1].get_position().y0, 
		0.025, 
		axes[0][-1].get_position().y1 - axes[1][-1].get_position().y0
	]) 
	plt.savefig("%s.pdf" % (STEM)) 
	plt.savefig("%s.png" % (STEM)) 

