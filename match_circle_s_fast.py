from __future__ import division
import numpy as np

def match_circle_s_fast(data1, data2, phase):

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
	
	T1_binn = np.zeros((2,len(T1_bin)))
	for i in range(len(T1_bin)):
		T1_binn[1,i] = np.sum(T1_bin[i])
		T1_binn[0,i] = T1_bin[i].size

	T2_binn = np.zeros((2,len(T2_bin)))
	for i in range(len(T2_bin)):
		T2_binn[1,i] = np.sum(T2_bin[i])
		T2_binn[0,i] = T2_bin[i].size

	T = np.zeros((2,len(T1_binn[0,:])))
	

	for i in range(len(T1_binn[1,:])):
		if T1_binn[0,i] != 0:
			T[0,i] = T1_binn[1,i]/T1_binn[0,i]
		else:
			T[0,i] = 0
		if T2_binn[0,i] != 0:
			T[1,i] = T2_binn[1,i]/T2_binn[0,i]
		else:
			T[1,i] = 0
		
	s_den = 0
	m = np.arange(0,360,1)
	T_m_a = np.fft.fft((T[0]))
	T_m_b = np.fft.fft((T[1]))
	T_m_a_sum = np.abs(np.sum(T_m_a))
	T_m_b_sum = np.abs(np.sum(T_m_b))
	lag = np.full(360, np.exp(-1j*phase))
	m_lag = np.multiply(m,lag)
	T_m_a = np.multiply(T_m_a, m)
	for m in range(360):
		s_den += m*((T_m_a_sum)**2+(T_m_b_sum)**2)

	if s_den != 0:
		s_m = (2*np.multiply(T_m_a,np.conj(T_m_b)))/s_den
		s = np.sum(np.multiply(s_m, m_lag))
	else:
		s=0
	
	return s
