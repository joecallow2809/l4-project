from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import random
from astropy.io import fits
from astrotools import healpytools as hpt
from circle_finder import circle_finder
from strip_auto import strip_finder
from load_file import load_file
from match_circle_r import match_circle_r


cmb_map = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")

NSIDE=2048
CMB_DIST = 14000
CELL_SIZE = 320

angle = random.uniform(0,np.pi/2)

strip_finder(cmb_map, circle_finder(CELL_SIZE, 0, 1), NSIDE)

circle_a = load_file('strip_a')
circle_b = load_file('strip_b')

x_corr_1 = match_circle(circle_a, circle_b)[1]

cmb_map_2 = hpt.rotate_map(cmb_map, np.array([0,0,1]), angle)

strip_finder(cmb_map_2, circle_finder(CELL_SIZE, 0, 1), NSIDE)

circle_a = load_file('strip_a')
circle_b = load_file('strip_b')

x_corr_2 = match_circle(circle_a, circle_b)[1]

cmb_map_3 = hpt.rotate_map(cmb_map_2, np.array([0,0,1]), -angle)

strip_finder(cmb_map_3, circle_finder(CELL_SIZE, 0, 1), NSIDE)

circle_a = load_file('strip_a')
circle_b = load_file('strip_b')

x_corr_3 = match_circle(circle_a, circle_b)[1]

fig1, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey = True, figsize = (30,16))

ax1.plot(np.arange(0,360,1), x_corr_1)
ax1.set_title('Auto-Correlation of Circle of CMB')
ax1.set_xlabel('Lag / $\circ$')
ax1.set_ylabel('X-Correlation')

ax2.plot(np.arange(0,360,1), x_corr_2)
ax2.set_title('Auto-Correlation of Rotated Circle of CMB')
ax2.set_xlabel('Lag / $\circ$')

ax3.plot(np.arange(0,360,1), x_corr_3)
ax3.set_title('Auto-Correlation of Rotated Circle of CMB, Rotated Back')
ax3.set_xlabel('Lag / $\circ$')

fig1.savefig('/opt/local/l4astro/rbbg94/figures/auto_corr_rot.png', overwrite = True)

fig2, (ax1, ax2, ax3) = plt.subplots(1,3, sharex=True, figsize=(30,8))

ax1.plot(np.arange(0,360,1), x_corr_1-x_corr_2)
ax1.set_title('Residuals between Normal and Rotated auto-correlations')
ax1.set_xlabel('Lag / $\circ$')
ax1.set_ylabel('Residuals')

ax2.plot(np.arange(0,360,1), x_corr_2-x_corr_3)
ax2.set_title('Residuals between Rotated and Rotated back auto-correlations')
ax2.set_xlabel('Lag / $\circ$')

ax3.plot(np.arange(0,360,1), x_corr_3-x_corr_1)
ax3.set_title('Residuals between Rotated back and Normal auto-correlations')
ax3.set_xlabel('Lag / $\circ$')

fig2.savefig('/opt/local/l4astro/rbbg94/figures/rot_residuals.png', overwrite = True)

fig3, (ax1, ax2, ax3) = plt.subplots(1,3,figsize = (20,15))

plt.axes(ax1)
hp.mollview(cmb_map, hold=True)
plt.axes(ax2)
hp.mollview(cmb_map_2, hold=True)
plt.axes(ax3)
hp.mollview(cmb_map_3, hold=True)

fig3.savefig('/opt/local/l4astro/rbbg94/figures/rotate.png', overwrite = True)

plt.show()
