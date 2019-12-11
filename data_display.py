from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('/opt/local/l4astro/rbbg94/downloads/torusres.dat')

data = data[0:50]

fig, ax = plt.subplots()

ax.plot(data[:,0], data[:,1])

plt.tight_layout()

plt.show()
