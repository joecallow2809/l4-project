from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

SMICA_MAP = hp.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")

NSIDE=2048

strip = hp.query_strip(NSIDE, np.radians(40), np.radians(50))

SMICA_MAP[strip] = SMICA_MAP.max()

hp.mollview(SMICA_MAP)

plt.show()
