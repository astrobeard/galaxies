""" 
This script prints the number of stars that migrate inward, outward, and at the 
same galactocentric radius as birth from the UWhydro simulation. 
""" 

import sys 
import os 
current = os.getcwd() 
os.chdir("../../") 
sys.path.append(os.getcwd()) 
os.chdir(current) 
from data import UWhydro 

def main(): 
	data = [list(i) for i in zip(
		UWhydro["rform"], 
		UWhydro["rfinal"] 
	)]
	inward = len(list(filter(lambda x: x[0] > x[1], data))) 
	outward = len(list(filter(lambda x: x[0] < x[1], data))) 
	same = len(list(filter(lambda x: x[0] == x[1], data))) 
	print("%d star particles migrated inward" % (inward)) 
	print("%d star particles migrated outward" % (outward)) 
	print("%d star particles at same galactocentric radius as birth" % (same)) 

if __name__ == "__main__": 
	main() 


