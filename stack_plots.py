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
from stack import stack

cmb_map_og = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")
NSIDE = hp.npix2nside(len(cmb_map_og))

apo = fits.open("/opt/local/l4astro/rbbg94/cmb_maps/mask_apo5.fits")

data = apo[1].data

mask_nest = data['GAL080'][:]

mask_ring = hp.pixelfunc.reorder(mask_nest, inp = 'nested', out = 'ring', n2r = True)

CMB_DIST = 14000
CELL_SIZE = 320

lon = np.array([165, 0, 0])
lat = np.array([90, 0, 0])

ang_rad = np.arange((1/360)*2*np.pi, np.pi/2, (2*np.pi)/360)

x_corr = np.zeros((3,len(ang_rad)), dtype=complex)

for plot in range(3):
	if lat[plot] == 90:
		cmb_map = np.multiply(cmb_map_og, mask_ring)
	else:
		cmb_map = rotate_to_top(np.multiply(cmb_map_og,mask_ring), lon[plot], lat[plot])

	for i in range(len(ang_rad)):
		strip_finder(cmb_map, ang_rad[i], NSIDE)
	
		circle_a = load_file('strip_a')
		circle_b = load_file('strip_b')
		x_corr[plot,i] = match_circle_s(circle_a, circle_b, 0, 720)
	if plot == 0:
		lat[plot+1] = lat[plot]-90
		lon[plot+1] = lon[plot]
		if lat[plot+1] < -90:
			lat[plot+1] = -(180+lat[plot+1])
			lon[plot+1] += 180
			if lon[plot+1] >= 360:
				lon[plot+1] -= 360
	if plot == 1:
		lon[plot+1] = lon[plot]+90
		if lon[plot+1] >= 360:
			lon[plot+1] -= 360

		
ang_rad = ang_rad*(360/(2*np.pi))

stack_corr = stack(x_corr[0], x_corr[1], x_corr[2])

fig, ((ax0,ax1),(ax2,ax3)) = plt.subplots(2,2, figsize=(24,16))

ax0.plot(ang_rad, x_corr[0])
ax0.set_title('Lon='+str(lon[0])+', Lat='+str(lat[0]))
ax0.set_xlabel('Angular Radius')
ax0.set_ylabel('X-Correlation: S')
ax0.legend(['Lag = 0$^\circ$'])
ax0.axhline(0, color = 'black')

ax1.plot(ang_rad, x_corr[1])
ax1.set_title('Lon='+str(lon[1])+', Lat='+str(lat[1]))
ax1.set_xlabel('Angular Radius')
ax1.set_ylabel('X-Correlation: S')
ax1.legend(['Lag = 0$^\circ$'])
ax1.axhline(0, color = 'black')

ax2.plot(ang_rad, x_corr[2])
ax2.set_title('Lon='+str(lon[2])+', Lat='+str(lat[2]))
ax2.set_xlabel('Angular Radius')
ax2.set_ylabel('X-Correlation: S')
ax2.legend(['Lag = 0$^\circ$'])
ax2.axhline(0, color = 'black')

ax3.plot(ang_rad, stack_corr)
ax3.set_title('Stacked')
ax3.set_xlabel('Angular Radius')
ax3.set_ylabel('X-Correlation: S')
ax3.legend(['Lag = 0$^\circ$'])
ax3.axhline(0, color = 'black')

fig.savefig('/opt/local/l4astro/rbbg94/figures/corr_0_stack_lon165_plane.png', overwrite = True)
plt.show()
