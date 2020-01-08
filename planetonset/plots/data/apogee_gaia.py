""" 
Subroutines for reading in and working with the APOGEE+Gaia data 
""" 

import numpy as np 
import os 

DIRS = os.path.abspath(__file__).split('/')[:-3] 
PATH = "/" 
for i in DIRS: 
	PATH += "%s/" % (i) 
GIANTS_FILE = "%sdata/APOGEE/dr16_giants.dat" % (PATH) 
DWARFS_FILE = "%sdata/APOGEE/dr16_dwarfs.dat" % (PATH) 

def giants(): 
	""" 
	Returns a dictionary containing the data on the dr16 halo giants. 
	""" 
	return _read_data_file(GIANTS_FILE) 

def dwarfs(): 
	""" 
	Returns a dictionary containing the data on the dr16 halo dwarfs. 
	""" 
	return _read_data_file(DWARFS_FILE) 

def whole(): 
	""" 
	Returns a dictionary containing the data on all dr16 halo stars. 
	""" 
	gi = giants() 
	dw = dwarfs() 
	data = {} 
	for i in gi.keys(): 
		data[i] = gi[i] + dw[i] 
	return data 

def _read_data_file(filename): 
	""" 
	Reads in one of the data files and returns a dictionary 
	""" 
	raw = np.genfromtxt(filename, delimiter = ',') 
	data = {} 
	data["m_h"] 			= [row[6] for row in raw] 
	data["m_h_err"] 		= [row[7] for row in raw] 
	data["alpha_m"] 		= [row[8] for row in raw] 
	data["alpha_m_err"] 	= [row[9] for row in raw] 
	return data 	

