from __future__ import division
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from astrotools import coord

def rotate_map(healpy_map, rotation_axis, rotation_angle):
    """
    Perform rotation of healpy map for given rotation axis and angle(s).

    :param healpy_map: healpix map to be rotated
    :param rotation_axis: rotation axis, either np.array([x, y, z]) or ndarray with shape (3, n)
    :param rotation_angle: rotation angle in radians, either float or array size n
    :return: rotated healpy map, same shape as input healpy map
    """
    nside = hp.get_nside(healpy_map)
    npix = hp.nside2npix(nside)
    _vecs = coord.rotate(hp.pix2vec(nside, np.arange(npix)), rotation_axis, -rotation_angle)
    _phi, _theta = hp.vec2ang(_vecs)
    return hp.get_interp_val(healpy_map, np.pi / 2. - _theta, _phi)



cmb_map = hp.fitsfunc.read_map("/opt/local/l4astro/rbbg94/cmb_maps/planck_data.fits")



NSIDE=2048

disc = hp.query_disc(NSIDE, (0,1,0), 0.1)

cmb_map[disc] = cmb_map.max()

cmb_map_2 = rotate_map(cmb_map, np.array([0,0,1]), np.pi/4)

fig, (ax1, ax2) = plt.subplots(1,2)


plt.sca(ax1)
hp.mollview(cmb_map, hold=True)
plt.sca(ax2)
hp.mollview(cmb_map_2, hold=True)

plt.show()
