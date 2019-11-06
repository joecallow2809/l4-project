from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from astrotools import healpytools as hpt

cmb_map = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")

NSIDE=2048


z_vec = np.array([0,0,1])

disc = hp.query_disc(NSIDE, z_vec, 5*(np.pi/180))

cmb_map[disc] = cmb_map.max()





hp.mollview(cmb_map, hold=True)


plt.show()
