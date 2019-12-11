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
from match_circle_s import match_circle_s
from rotate import rotate_to_top


cmb_map_og = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")
NSIDE = hp.npix2nside(len(cmb_map_og))

apo = fits.open("/opt/local/l4astro/rbbg94/cmb_maps/mask_apo5.fits")

data = apo[1].data

mask_nest = data['GAL080'][:]

mask_ring = hp.pixelfunc.reorder(mask_nest, inp = 'nested', out = 'ring', n2r = True)

CMB_DIST = 14000
CELL_SIZE = 320

lon = 207.8
lat = -56.3

cmb_map = np.multiply(cmb_map_og,mask_ring)
#cmb_map = rotate_to_top(np.multiply(cmb_map_og,mask_ring), lon, lat)

"""ang_rad = np.zeros(100)
x_corr = np.zeros(100)
a = -1

for i in range(len(ang_rad)):
	xi = (i+1)*1
	yi = (i+1)*0
	if np.sqrt(xi**2+yi**2) <= 87.5:
		strip_finder(cmb_map, circle_finder(CELL_SIZE, xi, yi), NSIDE)

		circle_a = load_file('strip_a')
		circle_b = load_file('strip_b')
	
		ang_rad[i] = (360/(2*np.pi))*circle_finder(CELL_SIZE, xi, yi)
		x_corr[i] = match_circle_r(circle_a, circle_b)[0]
		a += 1
	else:
		break

ang_rad = ang_rad[:a]
x_corr = x_corr[:a]

fig, ax = plt.subplots()

ax.plot(ang_rad, x_corr)
ax.set_title('Correlation of Circle of CMB: Cold Spot, domain direction')
ax.set_xlabel('Angular Radius')
ax.set_ylabel('X-Correlation')
ax.legend(['Lag = 0$^\circ$'])
ax.axhline(0, color = 'black')
plt.xticks(np.arange(0, 91, step=10))
plt.xlim(0,90)
plt.tight_layout()

fig.savefig('/opt/local/l4astro/rbbg94/figures/corr_pred_cold_spot__dom_dir_mask.png', overwrite = True)"""

"""strip_finder(cmb_map, 88.6*(2*np.pi/360), NSIDE)

circle_a = load_file('strip_a')
circle_b = load_file('strip_b')
	
x_corr = match_circle_r(circle_a, circle_b)

fig, ax = plt.subplots(figsize=(24,10))

ax.errorbar(np.arange(0,359,1), x_corr, yerr = 0, ecolor = 'red')
ax.set_title('Correlation of Circle of CMB: Cold-Spot')
ax.set_xlabel('Lag')
ax.set_ylabel('X-Correlation, r')
ax.legend([r'$\alpha$=87.9$^\circ$'])
ax.axhline(0, color = 'black')
plt.xticks(np.arange(0, 361, step=30))
plt.xlim(0,360)
plt.tight_layout()

fig.savefig('/opt/local/l4astro/rbbg94/figures/corr_cs_ang_88.6.png', overwrite = True)"""

ang_rad = np.arange((1/360)*2*np.pi, np.pi/2, (2*np.pi)/360)
x_corr = np.zeros(len(ang_rad), dtype=complex)
lag = 0
rad_lag = lag*np.pi/180

for i in range(len(ang_rad)):
	strip_finder(cmb_map, ang_rad[i], NSIDE)
	
	circle_a = load_file('strip_a')
	circle_b = load_file('strip_b')
	#x_corr[i] = match_circle_r(circle_a, circle_b)[0]
	x_corr[i] = match_circle_s(circle_a, circle_b, rad_lag, 720)
	
	
ang_rad = ang_rad*(360/(2*np.pi))


fig, ax = plt.subplots()

ax.plot(ang_rad, x_corr)
ax.set_title('Cross-Correlation of Circles of CMB: NGP_SGP')
ax.set_xlabel('Angular Radius')
ax.set_ylabel('X-Correlation: S ')
ax.legend(['Lag = '+str(lag)+'$^\circ$'])
ax.axhline(0, color = 'black')
plt.xticks(np.arange(0, 91, 10))
plt.xlim(0,90)
plt.tight_layout()

fig.savefig('/opt/local/l4astro/rbbg94/figures/corr_'+str(lag)+'ngp_sgp.png', overwrite = True)
plt.show()
