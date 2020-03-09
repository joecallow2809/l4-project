from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

ang_rad = np.linspace(1,89,89)
ratio = np.genfromtxt("/opt/local/l4astro/rbbg94/data/ngp_sim_err_2sims_0.csv", delimiter = ',')[:-1]/np.genfromtxt("/opt/local/l4astro/rbbg94/data/ngp_sim_err_2sims_1.csv", delimiter = ',')[:-1]

fig, ax = plt.subplots()
ax.plot(ang_rad, ratio)
plt.show()
