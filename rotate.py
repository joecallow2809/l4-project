from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from astrotools import healpytools as hpt


cmb_map = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")

NSIDE=2048

disc = hp.query_disc(NSIDE, (0,1,0), 0.1)

cmb_map[disc] = cmb_map.max()

cmb_map_2 = hpt.rotate_map(cmb_map, np.array([0,0,1]), np.pi/2)
cmb_map_3 = hpt.rotate_map(cmb_map_2, np.array([0,0,1]), -np.pi/2)

fig, (ax1, ax2, ax3) = plt.subplots(1,3)


plt.sca(ax1)
hp.mollview(cmb_map, hold=True)
plt.sca(ax2)
hp.mollview(cmb_map_2, hold=True)
plt.sca(ax3)
hp.mollview(cmb_map_3, hold=True)

fig.savefig('/opt/local/l4astro/rbbg94/figures/rotate.png', overwrite = True)

plt.show()
