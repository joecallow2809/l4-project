from __future__ import division
import numpy as np
import healpy as hp
import random
import matplotlib.pyplot as plt
from astropy.io import fits
from astrotools import healpytools as hpt
from strip import strip_finder
from load_file import load_file
from rotate import rotate_to_top


cmb_map = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")

NSIDE=hp.npix2nside(len(cmb_map))

apo = fits.open("/opt/local/l4astro/rbbg94/cmb_maps/mask_no_apo.fits")

data = apo[1].data

mask_nest = data['GAL090'][:]

mask_ring = hp.pixelfunc.reorder(mask_nest, inp = 'nested', out = 'ring', n2r = True)

cmb_map_masked = cmb_map*mask_ring

diff = cmb_map-cmb_map_masked

hp.mollview(cmb_map)
hp.mollview(cmb_map_masked)
print diff
hp.mollview(diff)
plt.show()
