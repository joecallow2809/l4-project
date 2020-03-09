from __future__ import division
from astropy.io import fits
import numpy as np


def load_file(fname, bins):
	hdul = fits.open('/opt/local/l4astro/rbbg94/cmb_maps/'+fname)
	circ_dat = hdul[1].data
	hdul.close()
	T = circ_dat['T']
	lon = circ_dat['long']
	lat = circ_dat['lat']
	index = circ_dat['index']
		
	circle = np.zeros((len(T),4))
	circle[:,3] = index
	circle[:,0] = T
	circle[:,1] = lon
	circle[:,2] = lat
	
	T_bin = np.zeros((bins, 2))
	for i in range(len(circle[:,1])):
		iy = int(circle[i,1])

		T_bin[iy,1] += circle[i,0]
		T_bin[iy,0] += 1
	
	zeroes = np.argwhere(T_bin[:,0]==0)
	np.put(T_bin[0], zeroes, 1)
	np.put(T_bin[1], zeroes, 0)
	T_binn = T_bin[:,1]/T_bin[:,0]

	return T_binn
