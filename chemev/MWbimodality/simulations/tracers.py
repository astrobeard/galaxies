""" 
Subroutines for setting up tracer particles in multizone simulations 
by tuning them to hydrodynamical simulation star particles. 
""" 

__all__ = ["UWhydro", "UWhydro_inward", "UWhydro_outward", "UWhydro_reverse"] 
import numpy as np 

def _get_bin_number(bins, val): 
	""" 
	Gets the bin number of a given value in a given binspace. 

	bins :: the array containing the bin edges. Assumed to be sorted. 
	val :: the number whose bin number needs to be found 

	Returns -1 if the value does not lie within the binspace 
	""" 
	for i in range(len(bins) - 1): 
		if bins[i] <= val <= bins[i + 1]: return i 
	return -1 

def _interpolate(x1, x2, y1, y2, x): 
	""" 
	Interpolate between two points (x1, y1) and (x2, y2) 

	Parameters 
	========== 
	x1 :: float 
		The x-coordinate of the first point 
	y1 :: float 
		The y-coordinate of the first point 
	x2 :: float 
		The x-coordinate of the second point 
	y2 :: float 
		The y-coordinate of the second point 
	x :: float 
		The x-coordinate of the point to interpolate to 

	Returns 
	======= 
	y :: float 
		The y-coordinate of the linearly extrapolated point 
	""" 
	try: 
		return (y2 - y1) / (x2 - x1) * (x - x1) + y1 
	except ZeroDivisionError: 
		return y2 


def _rand_range(start, stop): 
	r""" 
	Return a randomly-generated number in a given range. 

	Parameters 
	----------
	start : real number 
		The lower-bound of the range 
	stop : real number 
		The upper-bound of the range 

	Returns 
	-------
	x : real number 
		A pseudo-randomly generated number in the range [start, stop). 
	""" 
	return start + (stop - start) * np.random.random() 


class UWhydro(object): 

	""" 
	A callable object with tracer paticle data tuned to the UW hydro simulation 
	interpolating linearly between zone numbers. 
	""" 

	def __init__(self, time_bins, rad_bins, n_stars = 1, 
		filename = "tracers.out"): 
		self._time_bins = time_bins[:] 
		self._rad_bins = rad_bins[:] 
		self._n_stars = n_stars 
		self._zones, self._heights = self._analyze_radii() 
		self._file = open(filename, 'w') 
		self._file.write("# zone_origin\ttime_origin\tzone_final\tzfinal\n") 
		self.write = False 

	def __call__(self, zone, time, t, n = 0): 
		# tbin = _get_bin_number(self._time_bins, time) 
		# idx = np.random.randint(len(self._zones[zone][tbin])) 
		# final = self._zones[zone][tbin][idx] + np.random.random() 
		# init = zone + np.random.random() 
		# if self.write: 
		# 	self._file.write("%d\t%.3f\t%d\t%.3f\n" % (zone, time, final, 
		# 		self._heights[zone][tbin][idx])) 
		# else: 
		# 	pass 
		# def zones(t): 
		# 	if t < time: 
		# 		return 0 
		# 	elif t == time: 
		# 		return zone 
		# 	else: 
		# 		return int(_interpolate(time, self._time_bins[-1], init, final, 
		# 			t)) 
		# return zones 

		tbin = _get_bin_number(self._time_bins, time) 
		if t == time: 
			self._idx = np.random.randint(len(self._zones[zone][tbin])) 
			self._init = zone + np.random.random() 
			self._final = self._zones[zone][tbin][self._idx] + np.random.random() 
			if self.write: 
				self._file.write("%d\t%.2f\t%d\t%.3f\n" % (zone, time, 
					self._final, self._heights[zone][tbin][self._idx]))  
			else: 
				pass 
		else: 
			pass 
		if t < time: 
			return 0 
		elif t == time: 
			return zone 
		else: 
			return int(_interpolate(time, self._time_bins[-1], self._init, 
				self._final, t)) 


	def _analyze_radii(self): 
		print("Analyzing radii....") 
		from data import UWhydroparticles 
		zones = (len(self._rad_bins) - 1) * [None] 
		heights = (len(self._rad_bins) - 1) * [None] 
		for i in range(len(zones)): 
			zones[i] = len(self._time_bins) * [None] 
			heights[i] = len(self._time_bins) * [None] 
			for j in range(len(zones[i])): 
				zones[i][j] = [] 
				heights[i][j] = [] 
		for i in range(len(UWhydroparticles["tform"])): 
			tbin = _get_bin_number(self._time_bins, 
				UWhydroparticles["tform"][i]) 
			rbin = _get_bin_number(self._rad_bins, 
				UWhydroparticles["rform"][i]) 
			zones[rbin][tbin].append(_get_bin_number(self._rad_bins, 
				UWhydroparticles["rfinal"][i])) 
			heights[rbin][tbin].append(UWhydroparticles["zfinal"][i]) 
		for i in range(len(zones)): 
			for j in range(len(zones[i])): 
				if len(zones[i][j]) == 0: 
					# let it find something in a neighboring zone 
					if i > 0: 
						zones[i][j] += zones[i - 1][j][:] 
						heights[i][j] += heights[i - 1][j][:] 
					if i < len(zones) - 1: 
						zones[i][j] += zones[i + 1][j][:] 
						heights[i][j] += heights[i + 1][j][:] 
					if len(zones[i][j]) == 0: 
						zones[i][j].append(i) 
						heights[i][j].append(100) # ignore after the fact 
						# print("i = %d" % (i)) 
						# print(zones[i][j]) 
				else: continue 
		return [zones, heights] 


	def close_file(self): 
		self._file.close() 



class UWhydro_1event(UWhydro): 

	def __init__(self, time_bins, rad_bins, n_stars = 1, 
		filename = "tracers.out"): 
		super().__init__(time_bins, rad_bins, n_stars = n_stars, 
			filename = filename) 


	def __call__(self, zone, time, t, n = 0): 
		tbin = _get_bin_number(self._time_bins, time) 
		if t == time: 
			self._idx = np.random.randint(len(self._zones[zone][tbin])) 
			# self._init = zone 
			self._final = self._zones[zone][tbin][self._idx] 
			self._mig_time = _rand_range(time, 12.8) 
			if self.write: 
				self._file.write("%d\t%.2f\t%d\t%.3f\n" % (zone, time, 
					self._final, self._heights[zone][tbin][self._idx]))  
			else: 
				pass 
		else: 
			pass 
		if t == time: 
			return zone 
		elif t > self._mig_time: 
			return self._final 
		else: 
			return zone 
















# archived subclasses 

class UWhydro_zfilter(UWhydro): 

	def __init__(self, time_bins, rad_bins): 
		super().__init__(time_bins, rad_bins) 

	def _analyze_radii(self): 
		from data import UWhydroparticles_zfilter 
		zones = (len(self._rad_bins) - 1) * [None] 
		for i in range(len(zones)): 
			zones[i] = len(self._time_bins) * [None] 
			for j in range(len(zones[i])): 
				zones[i][j] = [] 
		for i in range(len(UWhydroparticles_zfilter["tform"])): 
			tbin = _get_bin_number(self._time_bins, 
				UWhydroparticles_zfilter["tform"][i]) 
			rbin = _get_bin_number(self._rad_bins, 
				UWhydroparticles_zfilter["rform"][i]) 
			zones[rbin][tbin].append(_get_bin_number(self._rad_bins, 
				UWhydroparticles_zfilter["rfinal"][i]))  
		for i in range(len(zones)): 
			for j in range(len(zones[i])): 
				if len(zones[i][j]) == 0: 
					zones[i][j].append(i) 
				else: continue 
		return zones 


class UWhydro_inward(UWhydro): 

	""" 
	A callable object with tracer particle data tuned to the UW hydro simulation 
	interpolating liearly between zone numbers, but selects only star particles 
	which migrate inward. 
	""" 

	def __init__(self, time_bins, rad_bins): 
		super().__init__(time_bins, rad_bins) 

	def __call__(self, zone, time): 
		tbin = _get_bin_number(self._time_bins, time) 
		possibilities = list(filter(lambda x: x <= zone, 
			self._zones[zone][tbin])) 
		if len(possibilities) > 0: 
			final = possibilities[np.random.randint(len(possibilities))] 
		else: 
			final = zone 
		final += np.random.random() 
		init = zone + np.random.random() 
		def zones(t): 
			if t < time: 
				return 0 
			elif t == time: 
				return zone 
			else: 
				return int(_interpolate(time, self._time_bins[-1], init, final, 
					t)) 
		return zones 
		

class UWhydro_outward(UWhydro): 

	""" 
	A callable object with tracer particle data tuned to the UW hydro simulation 
	interpolating liearly between zone numbers, but selects only star particles 
	which migrate outward. 
	""" 

	def __init__(self, time_bins, rad_bins): 
		super().__init__(time_bins, rad_bins) 

	def __call__(self, zone, time): 
		tbin = _get_bin_number(self._time_bins, time) 
		possibilities = list(filter(lambda x: x >= zone, 
			self._zones[zone][tbin])) 
		if len(possibilities) > 0: 
			final = possibilities[np.random.randint(len(possibilities))] 
		else: 
			final = zone 
		final += np.random.random() 
		init = zone + np.random.random() 
		def zones(t): 
			if t < time: 
				return 0 
			elif t == time: 
				return zone 
			else: 
				return int(_interpolate(time, self._time_bins[-1], init, final, 
					t)) 
		return zones 


class UWhydro_reverse(UWhydro): 

	""" 
	A callable object with tracer paticle data tuned to the UW hydro simulation 
	interpolating linearly between zone numbers, but with flipped initial and 
	final zone numbers 
	""" 

	def __init__(self, time_bins, rad_bins): 
		super().__init__(time_bins, rad_bins) 

	def _analyze_radii(self): 
		from data import UWhydroparticles 
		zones = (len(self._rad_bins) - 1) * [None] 
		for i in range(len(zones)): 
			zones[i] = len(self._time_bins) * [None] 
			for j in range(len(zones[i])): 
				zones[i][j] = [] 
		for i in range(len(UWhydroparticles["tform"])): 
			tbin = _get_bin_number(self._time_bins, 
				UWhydroparticles["tform"][i]) 
			rbin = _get_bin_number(self._rad_bins, 
				UWhydroparticles["rfinal"][i]) 
			zones[rbin][tbin].append(_get_bin_number(self._rad_bins, 
				UWhydroparticles["rform"][i]))  
		for i in range(len(zones)): 
			for j in range(len(zones[i])): 
				if len(zones[i][j]) == 0: 
					zones[i][j].append(i) 
				else: continue 
		return zones 

