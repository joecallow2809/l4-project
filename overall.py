from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from circle_finder import circle_finder
from strip import strip_finder
from load_file import load_file
from match_circle import match_circle

data = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")

NSIDE=2048
CMB_DIST = 14000
CELL_SIZE = 320

strip_finder(data, circle_finder(CELL_SIZE, 1, 0), NSIDE)

circle_a = load_file('strip_a')
circle_b = load_file('strip_b')

print match_circle(circle_a, circle_b)
