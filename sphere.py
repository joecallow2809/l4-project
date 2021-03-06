from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_aspect("equal")


theta = np.linspace(0, 2*np.pi, 360)
y = np.cos(theta)
x = np.sin(theta)

phi = 0
ax.plot(x, y*np.cos(phi), y*np.sin(phi)+0.03, color = 'b')
ax.plot(x, y*np.cos(phi), y*np.sin(phi)-0.03, color = 'b')

phi = np.pi/2
ax.plot(x, y*np.cos(phi)+0.03, y*np.sin(phi), color = 'r')
ax.plot(x, y*np.cos(phi)-0.03, y*np.sin(phi), color = 'r')

ax.view_init(elev=20, azim=7)

fig.savefig('/opt/local/l4astro/rbbg94/figures/circles_test.png', overwrite = True)
plt.show()
