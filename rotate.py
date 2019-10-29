from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits

class.hp.rotator.Rotator():

	data = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")

	NSIDE=2048

	disc = hp.query_disc(NSIDE, (0,1,0), 0.1)

	data[disc] = data.max()

	hp.mollview(data)

	plt.show()

	rotate_map_alms()

	hp.mollview(data)

	plt.show()

