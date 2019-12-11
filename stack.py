from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import inspect

def stack(*args):
	n = len(args)
	array = np.zeros(len(args[0]), dtype=complex)
	for arg in args:
		array += arg

	return array/n
