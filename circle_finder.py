from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

CMB_DIST = 14000 #Mpc

def circle_finder(domain_scale, x, y):
	d = domain_scale*(np.sqrt(x**2+y**2))
	ang_rad = np.arccos(d/(2*CMB_DIST))
	return ang_rad

domain_scale = 320
no_of_circles = int(CMB_DIST/domain_scale)
x=0
y=0

for i in range(no_of_circles):
	x+=1
	print np.degrees(circle_finder(domain_scale, x, y))



