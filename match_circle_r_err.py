from __future__ import division
import numpy as np

def match_circle_r(data1, data2):
	
	
	T1_bin = [[] for _ in range(360)]
	T2_bin = [[] for _ in range(360)]

	for i in range(len(data1[:,1])):
		iy = int(data1[i,1])
		if data1[i,0] != 0:
			T1_bin[iy].append(data1[i,0])
	T1_bin = np.array([np.array(T1i) for T1i in T1_bin])
	
			
	for i in range(len(data2[:,2])):
		iy = int(data2[i,1])
		if data2[i,0] != 0:      	  
			T2_bin[iy].append(data2[i,0])
	T2_bin = np.array([np.array(T2i) for T2i in T2_bin])
	
	T1_binn = np.zeros((len(T1_bin),2))
	for i in range(len(T1_bin)):
		T1_binn[i,1] = np.sum(T1_bin[i])
		T1_binn[i,0] = T1_bin[i].size

	T2_binn = np.zeros((len(T2_bin),2))
	for i in range(len(T2_bin)):
		T2_binn[i,1] = np.sum(T2_bin[i])
		T2_binn[i,0] = T2_bin[i].size

	T = np.zeros((2,len(T1_binn[:,0])))
	

	for i in range(len(T1_binn[:,1])):
		if T1_binn[i,0] != 0:
			T[0,i] = T1_binn[i,1]/T1_binn[i,0]
		else:
			T[0,i] = 0
		if T2_binn[i,0] != 0:
			T[1,i] = T2_binn[i,1]/T2_binn[i,0]
		else:
			T[1,i] = 0

	sum_val = np.zeros(360)
	ind_val = np.zeros((360, len(T[0])))
	circle1_tot = 0
	circle2_tot = 0

	for i in range(360):
		circle1_tot += T[0,i]**2
		circle2_tot += T[1,i]**2
		for j in range(len(T[0])):
		    if i + j >= 360:
		        x = (i + j - 360)
		    else:
		        x = (i + j)
		    ind_val[i][j] = T[0,x]*T[1,j]
		    sum_val[i] = sum_val[i] + ind_val[i][j]
		    
	if circle1_tot != 0 or circle2_tot != 0:
		norm_all = sum_val/(np.sqrt(circle1_tot)*np.sqrt(circle2_tot))    
	else:
		norm_all = np.arange(0,359,1)
	
	
	var = np.zeros((2, len(T1_bin)))
	for i in range(len(var[0])):
		if T1_bin[i].size != 0:
			var[0,i] = np.sum((T1_bin[i] - T[0,i])**2)/(T1_bin[i].size-1)
		else:
			var[0,i] = 0
		if T2_bin[i].size != 0:
			var[1,i] = np.sum((T2_bin[i] - T[1,i])**2)/(T2_bin[i].size-1)
		else:
			var[1,i] = 0
		
	std_dev = np.sqrt(var)
	std_err = np.zeros((2, len(std_dev[0])))
	for i in range(len(std_err)):
		if T1_bin[i].size != 0:
			std_err[0,i] = std_dev[0,i]/np.sqrt(T1_bin[i].size)
		else:
			std_err[0,i] = 0
		if T2_bin[i].size != 0:
			std_err[1,i] = std_dev[1,i]/np.sqrt(T2_bin[i].size)
		else:
			std_err[1,i] = 0

	sum_val_x = np.zeros(360)
	ind_val_x = np.zeros((360, len(T[0])))
	circle1_tot_x = 0
	circle2_tot_x = 0
	sum_val_y = np.zeros(360)
	ind_val_y = np.zeros((360, len(T[0])))
	circle1_tot_y = 0
	circle2_tot_y = 0

	for i in range(360):
		circle1_tot_x += (T[0,i]+std_err[0,i])**2
		circle2_tot_x += (T[1,i])**2
		circle1_tot_y += T[0,i]**2
		circle2_tot_y += (T[1,i]+std_err[1,i])**2
		for j in range(len(T[0])):
		    if i + j >= 360:
		        x = (i + j - 360)
		    else:
		        x = (i + j)
		    ind_val_x[i][j] = (T[0,x]+std_err[0,x])*(T[1,j])
		    sum_val_x[i] = sum_val_x[i] + ind_val_x[i][j]
		    ind_val_y[i][j] = T[0,x]*(T[1,j]+std_err[1,j])
		    sum_val_y[i] = sum_val_y[i] + ind_val_y[i][j]
		    
	if circle1_tot_x != 0 or circle2_tot_x != 0:
		x_err = sum_val_x/(np.sqrt(circle1_tot_x)*np.sqrt(circle2_tot_x))    
	else:
		x_err = np.zeros(len(norm_all))
	    
	if circle1_tot_y != 0 or circle2_tot_y != 0:
		y_err = sum_val_y/(np.sqrt(circle1_tot_y)*np.sqrt(circle2_tot_y))    
	else:
		y_err = np.zeros(len(norm_all))

	err = np.sqrt((x_err-norm_all)**2+(y_err-norm_all)**2)
	
	return norm_all, err
    
