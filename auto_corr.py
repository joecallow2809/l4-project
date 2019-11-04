from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from circle_finder import circle_finder
from strip_auto import strip_finder
from load_file import load_file
from match_circle import match_circle
from astrotools import healpytools as hpt

data = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")

NSIDE=2048
CMB_DIST = 14000
CELL_SIZE = 320

strip_finder(data, circle_finder(CELL_SIZE, 0, 0, 1), NSIDE)

circle_a = load_file('strip_a')
circle_b = load_file('strip_b')

x_corr = match_circle(circle_a, circle_b)[1]

fig, ax = plt.subplots()

ax.plot(np.arange(0,360,1), x_corr)
ax.set_title('Auto-Correlation of Circle of CMB')
ax.set_xlabel('Lag / $\circ$')
ax.set_ylabel('X-Correlation')

fig.savefig('/opt/local/l4astro/rbbg94/figures/auto_corr.png', overwrite = True)


plt.show()
