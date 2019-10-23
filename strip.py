from __future__ import division
import numpy as np
import healpy as hp
from circle_finder import circle_finder
import matplotlib.pyplot as plt


SMICA_MAP = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")

NSIDE=2048
CMB_DIST = 14000
CELL_SIZE = 320

ang_rad = circle_finder(CELL_SIZE, 1, 0)

ipix_strip_1 = hp.query_strip(NSIDE, ang_rad-(np.pi/360), ang_rad+(np.pi/360))
ipix_strip_2 = hp.query_strip(NSIDE, np.pi-ang_rad-(np.pi/360), np.pi-ang_rad+(np.pi/360))

SMICA_MAP[ipix_strip_1] = SMICA_MAP.max()
SMICA_MAP[ipix_strip_2] = SMICA_MAP.min()
hp.mollview(SMICA_MAP)

plt.show()


