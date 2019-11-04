from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from circle_finder import circle_finder
from strip import strip_finder
from load_file import load_file
from match_circle import match_circle
from astrotools import healpytools as hpt

data = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")

NSIDE = 2048
CMB_DIST = 14000
CELL_SIZE = 320

ang_rad = np.zeros(87)
x_corr = np.zeros(87)

for i in range(87):
	strip_finder(data, circle_finder(CELL_SIZE, 0, 0, i), NSIDE)

	circle_a = load_file('strip_a')
	circle_b = load_file('strip_b')
	
	ang_rad[i] = (360/(2*np.pi))*circle_finder(CELL_SIZE, 0, 0, i)
	x_corr[i] = match_circle(circle_a, circle_b)[1][29]

fig, ax = plt.subplots()

ax.plot(ang_rad, x_corr)
ax.set_title('Correlation of Circle of CMB')
ax.set_xlabel('Angular Radius')
ax.set_ylabel('X-Correlation')
ax.legend(['Lag = 30$^\circ$'])

fig.savefig('/opt/local/l4astro/rbbg94/figures/corr_30.png', overwrite = True)


plt.show()
