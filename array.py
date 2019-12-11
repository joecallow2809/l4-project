import numpy as np

array = [[] for _ in range(10)]
for n in range(4):
	for i in range(len(array)):
		array[i].append(i*n)

np_array = np.array([np.array(arrayi) for arrayi in array])

np_array = (np_array -1)**2

array = np.arange(0,360,1)
print len(array)
