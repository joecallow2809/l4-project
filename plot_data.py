from __future__ import division
import numpy as np
from matplotlib import pyplot as plt

ang_rad = np.linspace(1,89,89)
x_corr = np.genfromtxt('/opt/local/l4astro/rbbg94/data/ngp_corr_257.csv', dtype=complex, delimiter = ',')[:-1]
errors = np.genfromtxt('/opt/local/l4astro/rbbg94/data/ngp_lag_err.csv', dtype=complex, delimiter = ',')

fig, ax = plt.subplots(figsize = (14,10))
#ax.errorbar(ang_rad, x_corr, yerr = errors, color = 'k', ecolor = 'k', elinewidth = 0.8, capsize = 3)
ax.plot(ang_rad, x_corr, color = 'k')
plt.axhline(y=0, color = 'k')
plt.xticks(np.arange(0, 91, 10))
plt.xlim(0,90)
ax.set_xlabel(r'$\alpha/^\circ$')
ax.set_ylabel('$S$')
ax.annotate('Phase=0$^\circ$', xy = (0.8,0.95), xycoords = 'axes fraction', bbox=dict(facecolor='none', edgecolor='black'))
ax.fill_between(ang_rad, 2*errors, -2*errors, color = 'lightgrey')
ax.fill_between(ang_rad, errors, -errors, color = 'darkgrey')
plt.tight_layout()

#fig.savefig('/opt/local/l4astro/rbbg94/figures/corr_0_ngp_lag_err.png', overwrite = True)
plt.show()
