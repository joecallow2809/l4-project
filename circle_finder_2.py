from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

CMB_DIST = 14000
CELL_SIZE = 320

def draw_circle(x, y, r, coords, angles):
	for i in range(len(angles)):
		coords[0] = x + r*np.cos(angles)
		coords[1] = y + r*np.sin(angles)
	return coords 

coords_1 = np.zeros((2, 360))
coords_2 = np.zeros((2, 360))
angles = np.radians(np.arange(0,360, 1))

circle_1 = draw_circle(0, 0, CMB_DIST, coords_1, angles)
circle_2 = draw_circle(1*CELL_SIZE, 0, CMB_DIST, coords_2, angles)

crosses = np.zeros((2,2))

diff = float('inf')

for i in range(len(circle_1)):
	for n in range(len(circle_2)):
		x_diff = abs(circle_1[0,i] - circle_2
		

plt.plot(circle_1[0], circle_1[1])
plt.plot(circle_2[0], circle_2[1])
plt.plot(crosses[0,0], crosses[1,0], 'x')
plt.show()

