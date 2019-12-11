from __future__ import division
import numpy as np

def match_circle_s(data1, data2, phase, m_max):

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
	
	T1_b = np.zeros(len(T1), dtype=complex)
	T2_b = np.zeros(len(T2), dtype=complex)
	for i in range(len(T1)):
		T1_b[i] = np.exp(-1j*i)
		T2_b[i] = np.exp(-1j*i)
		
	
	s_num = 0
	s_den = 0
	for m in range(m_max):
		T1_b_m = T1_b**m
		T2_b_m = T2_b**m
		T1_m = T1*T1_b_m
		T2_m = T2*T2_b_m
		T2_m_star = np.conj(T2_m)
		s_num += (m/(2*np.pi**2))*np.sum(T1_m)*np.sum(T2_m_star)*np.exp(-1j*m*phase)
		s_den += (m/(4*np.pi**2))*((np.abs(np.sum(T1_m)))**2+(np.abs(np.sum(T2_m)))**2)
	
	if s_den != 0:
		s = s_num/s_den
	else:
		s=0
	return s
