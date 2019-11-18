from __future__ import division
import numpy as np
import healpy as hp
import random
import matplotlib.pyplot as plt
from astropy.io import fits
from astrotools import healpytools as hpt


cmb_map = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")

NSIDE=2048

apo = fits.open("/opt/local/l4astro/rbbg94/cmb_maps/mask_apo5.fits")

data = apo[1].data

mask = data['GAL090'][:]

print cmb_map, mask

cmb_map_masked = cmb_map*mask

print cmb_map_masked

hp.fitsfunc.write_map('/opt/local/l4astro/rbbg94/cmb_maps/my_mask.fits', cmb_map_masked, overwrite = True)
