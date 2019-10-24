from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from circle_finder import circle_finder
from strip import strip_finder
from load_file import load_file


NSIDE=2048
CMB_DIST = 14000
CELL_SIZE = 320

strip_finder(circle_finder(CELL_SIZE, 1, 0), CELL_SIZE, NSIDE)
circlea = load_file('circlea.fits')
circleb = load_file('circleb.fits')


