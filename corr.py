from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from astrotools import healpytools as hpt
from circle_finder import circle_finder
from strip import strip_finder
from load_file import load_file
from match_circle import match_circle
from rotate import rotate_to_top


cmb_map_og = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")

NSIDE = 2048
CMB_DIST = 14000
CELL_SIZE = 320

lon = 117.8
lat = 0

#cmb_map = cmb_map_og
cmb_map = rotate_to_top(cmb_map_og, lon, lat)

ang_rad = np.zeros(87)
x_corr = np.zeros(87)

for i in range(87):
	strip_finder(cmb_map, circle_finder(CELL_SIZE, 0, i), NSIDE)

	circle_a = load_file('strip_a')
	circle_b = load_file('strip_b')
	
	ang_rad[i] = (360/(2*np.pi))*circle_finder(CELL_SIZE, 0,i)
	x_corr[i] = match_circle(circle_a, circle_b)[1][270]

fig, ax = plt.subplots()

ax.plot(ang_rad, x_corr)
ax.set_title('Correlation of Circle of CMB: Galactic plane, Lon=117.8')
ax.set_xlabel('Angular Radius')
ax.set_ylabel('X-Correlation')
ax.legend(['Lag = 270$^\circ$'])

fig.savefig('/opt/local/l4astro/rbbg94/figures/corr_270_plane_117.png', overwrite = True)

"""strip_finder(cmb_map, circle_finder(CELL_SIZE, 3, 4), NSIDE)

circle_a = load_file('strip_a')
circle_b = load_file('strip_b')
	
x_corr = match_circle(circle_a, circle_b)[1]

fig, ax = plt.subplots()

ax.plot(np.arange(0,360,1), x_corr)
ax.set_title('Correlation of Circle of CMB')
ax.set_xlabel('Lag')
ax.set_ylabel('X-Correlation')"""



"""ang_rad = np.arange(0, (2*np.pi/360)*90, (2*np.pi)/360)
x_corr = np.zeros(len(ang_rad))

for i in range(len(ang_rad)):
	strip_finder(cmb_map, ang_rad[i], NSIDE)
	
	circle_a = load_file('strip_a')
	circle_b = load_file('strip_b')
	x_corr[i] = match_circle(circle_a, circle_b)[1][0]
	
ang_rad = ang_rad*(360/(2*np.pi))	

fig, ax = plt.subplots()

ax.plot(ang_rad, x_corr)
ax.set_title('Correlation of Circle of CMB')
ax.set_xlabel('Angular Radius')
ax.set_ylabel('X-Correlation')"""

plt.show()
