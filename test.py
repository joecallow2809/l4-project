from __future__ import division
import numpy as np
import time



array = np.full((100,50),2)
array_2 = np.arange(0,10,0.1)

start = time.time()
array = np.swapaxes(array, 0,1)
array_3 = array**array_2
array_3 = np.swapaxes(array_3, 0,1)
print time.time()-start
print array_3

array = np.full((100,50),2)
start = time.time()
for i in range(5):
		array[:,i] = array[:,i]**array_2
print time.time()-start
print array

