
import matplotlib.pyplot as plt 
from matplotlib import colors, cm 
import plots 
plots.mpltoolkit.load_mpl_presets() 
from astropy.io import fits 
import numpy as np 
import math as m 
import vice 
from vice.yields.presets import JW20 
vice.yields.sneia.settings['fe'] *= 10**0.1 
import sys 
sys.path.append("/Users/astrobeard/Work/Research/VICErepos/VICE/migration") 
import src 

OUTPUTSDIR = "/Users/astrobeard/Work/Research/VICErepos/VICE/migration/outputs" 
FULL = "%s/high-resolution/2Gyr/diffusion/insideout" % (OUTPUTSDIR) 
SIMPLE = "%s/high-resolution/2Gyr/post-process/insideout" % (OUTPUTSDIR) 
STEM = "yar_insideout_highres" 
# FULL = "../../simulations/paper_noburst" 
# SIMPLE = "../../simulations/paper_noburst_simple"  
ZONE_WIDTH = 0.1 
CMAP = "jet" 
ZONE_MIN = int(7 / ZONE_WIDTH) 
ZONE_MAX = int((9 - ZONE_WIDTH) / ZONE_WIDTH) 

# cm1 = colors.LinearSegmentedColormap.from_list("MyCmap", ["r", "b"]) 
# cnorm = colors.Normalize(vmin = 0, vmax = 15) 
# cpick = cm.ScalarMappable(norm = cnorm, cmap = cm1) 
# cpick.set_array([]) 



def setup_axes(): 
	fig = plt.figure(figsize = (15, 6)) 
	axes = 3 * [None] 
	for i in range(3): 
		axes[i] = fig.add_subplot(131 + i, facecolor = "white") 
	plt.subplots_adjust(top = 0.8, bottom = 0.16, right = 0.98, left = 0.09) 
	pos0 = axes[0].get_position() 
	pos1 = axes[1].get_position() 
	pos1.x1 -= pos1.x0 - pos0.x1 
	pos1.x0 = pos0.x1 
	axes[1].set_position(pos1) 
	axes[0].set_ylabel("[O/Fe]") 
	axes[2].set_ylabel(r"$\propto\dot{N}_\text{Ia}$ [Gyr$^{-1}$]") 
	axes[2].set_xlabel("Time [Gyr]") 
	plt.setp(axes[1].get_yticklabels(), visible = False) 

	# use dummy axes to put x-axis label between the left and middle panels 
	dummy = fig.add_subplot(121, facecolor = "white", zorder = -1) 
	posd = dummy.get_position() 
	posd.x0 = pos0.x0 
	posd.x1 = pos1.x1 
	dummy.set_position(posd) 
	dummy.set_xlabel("Age [Gyr]", labelpad = 30) 
	plt.setp(dummy.get_xticklabels(), visible = False) 
	plt.setp(dummy.get_yticklabels(), visible = False) 

	for i in range(2): 
		axes[i].set_xlim([-1, 14]) 
		axes[i].set_ylim([-0.2, 0.5]) 
	axes[2].set_xlim([-1, 14]) 
	axes[2].set_ylim([-0.1, 2.1])
	axes[0].text(1, 0.35, "Final Timestep", fontsize = 20) 
	axes[1].text(1, 0.35, r"$\propto$ Age", fontsize = 20) 
	# axes[0].set_title("simple = True", fontsize = 20) 
	# axes[1].set_title("simple = False", fontsize = 20) 
	axes.append(dummy) 

	return axes 


def plot_relation(ax, output): 
	cmap = plt.get_cmap(CMAP) 
	stars = output.stars.filter("zone_final", ">=", ZONE_MIN) 
	stars = stars.filter("zone_final", "<=", ZONE_MAX) 
	stars = stars.filter("zfinal", ">=", -3.) 
	stars = stars.filter("zfinal", "<=", 3.) 
	stars = stars.filter("mass", ">=", 1.) 
	colors = [ZONE_WIDTH * (i + 0.5) for i in stars["zone_origin"]] 
	# colors = [cpick.to_rgba(ZONE_WIDTH * (i + 0.5)) 
		# for i in stars["zone_origin"]]
	return ax.scatter(stars["age"], stars["[O/Fe]"], c = colors, s = 0.1, 
		cmap = cmap, vmin = 0, vmax = 15) 
		# vmin = 0, vmax = 15)  


def zheights(name): 
	if "VICE" in name: 
		raw = np.genfromtxt("%s_analogdata.out" % (name)) 
	else: 
		raw = np.genfromtxt("%s_extra_tracer_data.out" % (name)) 
	return [row[-1] for row in raw] 


def feuillet2018_data(ax): 
	raw = np.genfromtxt("age_alpha.dat") 
	ofe = len(raw) * [0.] 
	age = len(raw) * [0.] 
	ofe_disp = len(raw) * [0.] 
	age_disp = len(raw) * [None] 
	age_disp = [len(raw) * [0.], len(raw) * [0.]] 
	for i in range(len(raw)): 
		ofe[i] = (raw[i][0] + raw[i][1]) / 2 
		ofe_disp[i] = (raw[i][1] - raw[i][0]) / 2 
		age[i] = 10**(raw[i][2] - 9) # -9 converts from log(age) in yr to Gyr 
		age_disp[0][i] = age[i] - 10**(raw[i][2] - raw[i][3] - 9) 
		age_disp[1][i] = 10**(raw[i][2] + raw[i][3] - 9) - age[i] 
	ax.errorbar(age, ofe, xerr = age_disp, yerr = ofe_disp, 
		c = plots.mpltoolkit.named_colors()["black"], linestyle = "None") 


def feuillet2019_data(ax): 
	raw = fits.open(
		"Feuillet2019_alpha_ages/ELEM_GAUSS_AGE_07_09_00_05_alpha.fits") 
	ofe = len(raw[1].data) * [0.] 
	age = len(raw[1].data) * [0.] 
	ofe_disp = len(raw[1].data) * [0.] 
	age_disp = [len(raw[1].data) * [0.], len(raw[1].data) * [0.]] 
	for i in range(len(raw[1].data)): 
		if raw[1].data['nstars'][i] > 15: 
			ofe[i] = (raw[1].data['bin_ab'][i] + 
				raw[1].data['bin_ab_max'][i]) / 2 
			ofe_disp[i] = (raw[1].data['bin_ab_max'][i] - 
				raw[1].data['bin_ab'][i]) / 2 
			age[i] = 10**(raw[1].data['mean_age'][i] - 9) 
			# age_disp[i] = 10**(raw[1].data['age_disp'][i] - 9) 
			age_disp[0][i] = age[i] - 10**(raw[1].data['mean_age'][i] - 
				raw[1].data['age_disp'][i] - 9) 
			age_disp[1][i] = 10**(raw[1].data['mean_age'][i] + 
				raw[1].data['age_disp'][i] - 9) - age[i] 
		else: 
			ofe[i] = float("nan") 
			age[i] = float("nan") 
			ofe_disp[i] = float("nan") 
			age_disp[0][i] = float("nan") 
			age_disp[1][i] = float("nan") 
	ax.errorbar(age, ofe, xerr = age_disp, yerr = ofe_disp, 
		c = plots.mpltoolkit.named_colors()["black"], linestyle = "None") 


def ia_rate_proxies(zone, prefactor = 1): 
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
				vice.yields.ccsne.settings['fe'] * zone.history['sfr'][i] * 1e9 
			) 
			proxies[i] += zone.history["z(fe)"][i] * zone.history['sfr'][i] * (
				1 + mir.eta - 0.4) * 1e9 
			proxies[i] /= zone.history["mass(fe)"][i] 
			proxies[i] *= prefactor 
			if proxies[i] < 0: proxies[i] = 0 
		else: pass 
	return proxies 


def plot_ia_rate_proxies(ax, output, linestyle = '-', label = True): 
	radii = [15, 10, 5] 
	colors = ["blue", "red", "black"]  
	# colors = ["black", "red", "blue"] 
	prefactors = [1.3, 1, 1] 
	for i in range(len(radii)): 
		kwargs = {
			"c": 			colors[i], 
			"linestyle": 	linestyle 
		} 
		if label: kwargs["label"] = "%g kpc" % (radii[i]) 
		zone = int(radii[i] / ZONE_WIDTH) 
		ax.plot(output.zones["zone%d" % (zone)].history["time"][:-1], 
			ia_rate_proxies(output.zones["zone%d" % (zone)], 
				prefactor = prefactors[i]), 
			**kwargs) 
	if label: 
		leg = ax.legend(loc = plots.mpltoolkit.mpl_loc("upper right"), 
			ncol = 1, frameon = False, bbox_to_anchor = (0.99, 0.99), 
			handlelength = 0, fontsize = 22)  
		for i in range(len(leg.legendHandles)): 
			leg.get_texts()[i].set_color(colors[i]) 
			leg.legendHandles[i].set_visible(False) 
	else: pass 


if __name__ == "__main__": 
	axes = setup_axes() 
	full = vice.output(FULL) 
	simple = vice.output(SIMPLE) 
	full.stars["zfinal"] = zheights(FULL)[:full.stars.size[0]] 
	simple.stars["zfinal"] = zheights(SIMPLE)[:simple.stars.size[0]] 
	sc = plot_relation(axes[0], simple) 
	plot_relation(axes[1], full) 
	plot_ia_rate_proxies(axes[2], full) 
	plot_ia_rate_proxies(axes[2], simple, linestyle = ':', label = False) 
	# for i in range(2): feuillet2018_data(axes[i]) 
	for i in range(2): feuillet2019_data(axes[i]) 
	# cbar_ax = plots.mpltoolkit.append_axes(axes[-1], loc = "top", 
	# 	size = "8%") 
	cbar_ax = plt.gcf().add_axes([0.1, 0.8, 0.5, 0.05]) 
	cbar = plt.colorbar(sc, cax = cbar_ax, pad = 0.0, 
		orientation = "horizontal")  
	cbar.set_label(r"$R_\text{gal}$ of birth [kpc]", labelpad = 10) 
	cbar.set_ticks(range(2, 16, 2))  
	cbar_ax.xaxis.set_ticks_position("top") 
	cbar_ax.xaxis.set_label_position("top") 
	cbar_ax.set_position([ 
		axes[0].get_position().x0, 
		axes[0].get_position().y1, 
		axes[1].get_position().x1 - axes[0].get_position().x0, 
		0.05 
	]) 
	# stem = "young_alpha_rich" 
	plt.savefig("%s.pdf" % (STEM)) 
	plt.savefig("%s.png" % (STEM)) 
	# stem = "young_alpha_rich_%s" % (CMAP) 
	# plt.savefig("%s.pdf" % (stem)) 
	# plt.savefig("%s.png" % (stem)) 


