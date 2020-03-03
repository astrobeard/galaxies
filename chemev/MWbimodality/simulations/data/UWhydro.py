""" 
This file reads in the UWhydro data and places it in a dictionary of 
formation and final radii along with formation times. 
""" 

import numpy as np 
import vice 
import os 


DIRS = os.path.abspath(__file__).split('/')[:-3] 
PATH = "/" 
for i in DIRS: 
	PATH += "%s/" % (i) 
FILE = "%sdata/UWhydro_modded.dat" % (PATH) 
RAW = np.genfromtxt(FILE)  
COLS = [1, 2, 4, 5, 6, 7, 8] 
LABELS = ["tform", "rform", "rfinal", "zfinal", "v_r", "v_phi", "v_z"] 
COPY = {} 
for i in range(len(COLS)): 
	COPY[LABELS[i]] = [row[COLS[i]] for row in RAW] 
UWhydro = vice.dataframe(COPY) 

FLTRD = list(filter(lambda x: abs(x[5]) <= 3 and abs(x[8]) <= 50, RAW)) 
COPY = {} 
for i in range(len(COLS)): 
	COPY[LABELS[i]] = [row[COLS[i]] for row in FLTRD]
UWhydro_zfilter = vice.dataframe(COPY) 

