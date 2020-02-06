
import numpy as np 
import vice 
import os 

raw = np.genfromtxt("%s/UWhydro_particles.dat" % (
	os.path.dirname(os.path.abspath(__file__)))) 
cols = [1, 2, 4, 5, 6, 7, 8] 
labels = ["tform", "rform", "rfinal", "zfinal", "v_r", "v_phi", "v_z"] 
copy = {} 
for i in range(len(cols)): 
	copy[labels[i]] = [row[cols[i]] for row in raw] 
UWhydro = vice.dataframe(copy) 

fltrd = list(filter(lambda x: abs(x[5]) <= 3 and abs(x[8]) <= 50, raw)) 
copy = {} 
for i in range(len(cols)): 
	copy[labels[i]] = [row[cols[i]] for row in fltrd] 
UWhydro_zfilter = vice.dataframe(copy) 

