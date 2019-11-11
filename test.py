from __future__ import division
import numpy as np
import healpy as hp
import random
import matplotlib.pyplot as plt
from astropy.io import fits
from astrotools import healpytools as hpt


cmb_map = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")

NSIDE=2048


y_vec = np.array([0,1,0])

disc = hp.query_disc(NSIDE, y_vec, 5*(np.pi/180))

cmb_map[disc] = cmb_map.max()

angle = random.uniform(0, np,pi/2)

cmb_map_2 = hpt.rotate_map(cmb_map, np.array([0,0,1]), angle)
cmb_map_3 = hpt.rotate_map(cmb_map_2, np.array([0,0,1]), -angle)

fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize = (20,15))

plt.axes(ax1)
hp.mollview(cmb_map, min=-0.000815966, hold=True)
plt.axes(ax2)
hp.mollview(cmb_map_2, min=-0.000815966, hold=True)
plt.axes(ax3)
hp.mollview(cmb_map_3, min=-0.000815966, hold=True)

fig.savefig('/opt/local/l4astro/rbbg94/figures/rotate.png', overwrite = True)

plt.show()
