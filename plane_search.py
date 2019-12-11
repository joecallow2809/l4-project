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
RANGE = 9

lon = 60
lat = 0

ang_rad = np.arange((1/360)*2*np.pi, np.pi/2, (2*np.pi)/360)
x_corr = np.zeros(len(ang_rad), dtype=complex)
x_corrs = np.zeros(shape=(360, RANGE))

for i in range(RANGE):
	
	cmb_map = rotate_to_top(np.multiply(cmb_map_og,mask_ring), lon, lat)

	for i in range(len(ang_rad)):
		strip_finder(cmb_map, ang_rad[i], NSIDE)
	
		circle_a = load_file('strip_a')
		circle_b = load_file('strip_b')
		x_corr[i] = match_circle_s(circle_a, circle_b, 0, 720)
	
	ang_rad_deg = ang_rad*(360/(2*np.pi))	

	fig, ax = plt.subplots()

	ax.plot(ang_rad_deg, x_corr)
	ax.set_title('Correlation of Circle of CMB: Galactic Plane, lon='+str(lon))
	ax.set_xlabel('Angular Radius')
	ax.set_ylabel('X-Correlation: S')
	ax.legend(['Lag = 0$^\circ$'])
	ax.axhline(0, color = 'black')
	plt.xticks(np.arange(0, 91, 10))
	plt.xlim(0,90)
	plt.tight_layout()

	fig.savefig('/opt/local/l4astro/rbbg94/figures/corr_0_plane_'+str(lon)+'.png', overwrite = True)

	lon += 180/12


