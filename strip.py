from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

SMICA_MAP = hp.read_map("COM_CMB_IQU-smica_2048_R3.00_full.fits")

NSIDE=2048

ipix_strip = hp.query_strip(NSIDE, np.radians(40), np.radians(50))

SMICA_MAP[ipix_strip] = SMICA_MAP.max()
hp.mollview(SMICA_MAP)

plt.show()