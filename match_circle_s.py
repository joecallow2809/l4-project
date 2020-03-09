from __future__ import division
import numpy as np
import time

def match_circle_s(data1, data2, phase, bins, m_max):

	T_m_1 = np.full((m_max,bins), np.exp(-1j), dtype = complex)
	T_m_2 = np.full((m_max,bins), np.exp(-1j), dtype = complex)
	m = np.arange(0,m_max,1)
	phi = np.arange(0,2*np.pi,2*np.pi/bins)
	T_m_1 = T_m_1**phi
	T_m_2 = T_m_2**phi
	for i in range(bins):
		T_m_1[:,i] = T_m_1[:,i]**m
		T_m_2[:,i] = T_m_2[:,i]**m

	T_m_1 = (1/(2*np.pi))*np.sum(data1*T_m_1, axis = 1)
	T_m_2 = (1/(2*np.pi))*np.sum(data2*T_m_2, axis = 1)
	T_m_2_star = np.conj(T_m_2)
	phase_factor = (np.full(m_max, np.exp(-1j*phase), dtype = complex))**m
	s_num = 2*np.sum(m*T_m_1*T_m_2_star*phase_factor)
	s_den = np.sum(m*((np.abs(T_m_1))**2+(np.abs(T_m_2))**2))

	if s_den != 0:
		s = s_num/s_den
	else:
		s=0
	
	return s
