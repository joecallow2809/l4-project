from __future__ import division
import numpy as np

def match_circle_r(data1, data2):

	T1_bin = np.zeros((360, 2))
	T2_bin = np.zeros((360, 2))
	for i in range(len(data1[:,1])):
		iy = int(data1[i,1])
		if data1[i,0] != 0:
			T1_bin[iy,1] += data1[i,0]
			T1_bin[iy,0] += 1
		
	for i in range(len(data2[:,2])):
		iy = int(data2[i,1])
		if data2[i,0] != 0:      	  
			T2_bin[iy,1] += data2[i,0]
			T2_bin[iy,0] += 1

	T1 = np.zeros(len(T1_bin[:,0]))
	T2 = np.zeros(len(T2_bin[:,0]))
	for i in range(len(T1_bin[:,1])):
		if T1_bin[i,0] != 0:
			T1[i] = T1_bin[i,1]/T1_bin[i,0]
		else:
			T1[i] = 0

	for i in range(len(T2_bin[:,1])):
		if T2_bin[i,0] != 0:
			T2[i] = T2_bin[i,1]/T2_bin[i,0]
		else:
			T2[i] = 0

	sum_val = np.zeros(360)
	ind_val = np.zeros((360, len(T1)))
	circle1_tot = 0
	circle2_tot = 0

	for i in range(360):
		circle1_tot += T1[i]**2
		circle2_tot += T2[i]**2
		for j in range(len(T1)):
		    if i + j >= 360:
		        x = (i + j - 360)
		    else:
		        x = (i + j)
		    ind_val[i][j] = T1[x]*T2[j]
		    sum_val[i] = sum_val[i] + ind_val[i][j]
		    
	if circle1_tot != 0 or circle2_tot != 0:
		norm_all = sum_val/(np.sqrt(circle1_tot)*np.sqrt(circle2_tot))    
	else:
		norm_all = np.arange(0,359,1)

	return norm_all
    
