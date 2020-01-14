""" 
This file reads in the UWhydro data and places it in a dictionary of 
formation and final radii along with formation times. 
""" 

import numpy as np 
import sys 
import os 

DIRS = os.path.abspath(__file__).split('/')[:-3] 
PATH = "/" 
for i in DIRS: 
	PATH += "%s/" % (i) 
FILE = "%sdata/UWhydro_particles.dat" % (PATH) 
RAW = np.genfromtxt(FILE) 
DATA = {} 
DATA["tform"] 		= [row[1] for row in RAW] 
DATA["rform"] 		= [row[2] for row in RAW] 
DATA["rfinal"] 		= [row[4] for row in RAW] 

