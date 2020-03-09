from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

SMICA_MAP = hp.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")

NSIDE=2048

strip_1 = hp.query_strip(NSIDE, np.radians(40), np.radians(41))

strip_2 = hp.query_strip(NSIDE, np.radians(130), np.radians(131))

SMICA_MAP[strip_1] = SMICA_MAP.max()
SMICA_MAP[strip_2] = SMICA_MAP.max()

hp.mollview(SMICA_MAP, title = '', cbar = False)

plt.savefig("/opt/local/l4astro/rbbg94/figures/cmb_map.png")

plt.show()
