from __future__ import division
from astropy.io import fits
import numpy as np


def load_file(fname):
    hdul = fits.open(fname)
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

    
    return circle
