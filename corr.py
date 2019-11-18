from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from astrotools import healpytools as hpt
from circle_finder import circle_finder
from strip import strip_finder
from load_file import load_file
from match_circle_r import match_circle_r
from rotate import rotate_to_top


cmb_map_og = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")

apo = fits.open("/opt/local/l4astro/rbbg94/cmb_maps/mask_apo5.fits")

data = apo[1].data

mask = data['GAL090'][:]

NSIDE = 2048
CMB_DIST = 14000
CELL_SIZE = 320

lon = 352
lat = -54

#cmb_map = cmb_map_og
cmb_map = rotate_to_top(cmb_map_og, lon, lat)

ang_rad = np.zeros(100)
x_corr = np.zeros(100)
a = -1

for i in range(len(ang_rad)):
	xi = (i+1)*2
	yi = (i+1)*3
	if np.sqrt(xi**2+yi**2) <= 87.5:
		strip_finder(cmb_map, circle_finder(CELL_SIZE, xi, yi), NSIDE)

		circle_a = load_file('strip_a')
		circle_b = load_file('strip_b')
	
		ang_rad[i] = (360/(2*np.pi))*circle_finder(CELL_SIZE, xi, yi)
		x_corr[i] = match_circle_r(circle_a, circle_b)[1][0]
		a += 1
	else:
		break

ang_rad = ang_rad[:a]
x_corr = x_corr[:a]

fig, ax = plt.subplots()

ax.plot(ang_rad, x_corr)
ax.set_title('Correlation of Circle of CMB: Cold Spot, NGP domain direction')
ax.set_xlabel('Angular Radius')
ax.set_ylabel('X-Correlation')
ax.legend(['Lag = 0$^\circ$'])
ax.axhline(0, color = 'black')

fig.savefig('/opt/local/l4astro/rbbg94/figures/corr_pred_cold_spot.png', overwrite = True)

"""strip_finder(cmb_map, 87.9, NSIDE)

circle_a = load_file('strip_a')
circle_b = load_file('strip_b')
	
x_corr = match_circle_r(circle_a, circle_b)[1]

fig, ax = plt.subplots()

ax.plot(np.arange(0,360,1), x_corr)
ax.set_title('Correlation of Circle of CMB: NGP-SGP')
ax.set_xlabel('Lag')
ax.set_ylabel('X-Correlation')
ax.legend([r'$\alpha$=87.9$^\circ$'])
ax.axhline(0, color = 'black')

fig.savefig('/opt/local/l4astro/rbbg94/figures/corr_ngp_ang_87.9.png', overwrite = True)"""

"""ang_rad = np.arange(0, np.pi/2, (2*np.pi)/360)
x_corr = np.zeros(len(ang_rad))

for i in range(len(ang_rad)):
	strip_finder(cmb_map, ang_rad[i], NSIDE)
	
	circle_a = load_file('strip_a')
	circle_b = load_file('strip_b')
	x_corr[i] = match_circle_r(circle_a, circle_b)[1][0]
	
ang_rad = ang_rad*(360/(2*np.pi))	

fig, ax = plt.subplots()

ax.plot(ang_rad, x_corr)
ax.set_title('Correlation of Circle of CMB: Lon=352$^\circ$, Lat=-54$^\circ$')
ax.set_xlabel('Angular Radius')
ax.set_ylabel('X-Correlation')
ax.legend(['Lag = 0$^\circ$'])
ax.axhline(0, color = 'black')

fig.savefig('/opt/local/l4astro/rbbg94/figures/corr_0_lon352_latminus54.png', overwrite = True)"""
plt.show()
