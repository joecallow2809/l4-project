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
import time

start = time.time()

CMB_DIST = 14000
CELL_SIZE = 320

lon = 207.8
lat = -56.3

no_sims = 100
bins = 360
error_phases = 360

ang_rad = np.linspace((1/360)*2*np.pi, np.pi/2, 90)
no_circles = len(ang_rad)
x_corr = np.zeros((no_sims,no_circles,error_phases), dtype=complex)

for i in xrange(no_sims):
	cmb_map = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/sims/dx12_v3_smica_cmb_mc_000"+str(i).zfill(2)+"_raw.fits")
	NSIDE = hp.npix2nside(len(cmb_map))

	for j in xrange(no_circles):
		rad_lag = 0
		strip_finder(cmb_map, ang_rad[j], NSIDE)
		circle_a = load_file('strip_a', bins)
		circle_b = load_file('strip_b', bins)
		for k in xrange(error_phases):
			x_corr[i,j,k] = match_circle_s(circle_a, circle_b, rad_lag, bins, 720)
			rad_lag += 2*np.pi/bins
			print k, time.time() - start

mean = np.zeros(no_circles)
std_dev = np.zeros(no_circles)

for i in xrange(error_phases):
	for j in xrange(no_circles):
		sum_sq = 0
		mean[j] = np.mean(x_corr[:,j,i])
		for k in xrange(no_sims):
			sum_sq += (x_corr[k,j,i]-mean[j])**2
		std_dev[j] = np.sqrt(sum_sq/(no_sims-1))
	np.savetxt("/opt/local/l4astro/rbbg94/data/ngp_sim_err_"+str(no_sims)+"_sims_"+str(i)+".csv", std_dev, delimiter = ',')
