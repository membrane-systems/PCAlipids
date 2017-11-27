import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Pool
import time


def estimated_autocorrelation(x):
	n = len(x)
	variance = x.var()
	x = x-x.mean()
	r = np.correlate(x, x, mode = 'full')[-n:]
	assert np.allclose(r, np.array([(x[:n-k]*x[-(n-k):]).sum() for k in range(n)]))
	result = r/(variance*(np.arange(n, 0, -1)))
	return result


def main(filename):
	file = open(filename, 'r')
	data = []
	for line in file:
	    if line.find('@') == -1 and line.find('&') == -1:
	        data.append(np.float(line.split()[0]))
	    if line.find('&') != -1:
	        break
	file.close()
	data = np.array(data, dtype = np.float64)
	data_1 = data[:len(data) // 128]
	R = estimated_autocorrelation(data_1)
	for i in range(1, 128):
	    data_1 = data[len(data) * i // 128:len(data) * (i + 1) // 128]
	    R += estimated_autocorrelation(data_1)
	R /= 128
	T = []
	for i in range(1, len(R) + 1):
	    T.append(0.01 * i)
	print(filename)
	return T, R

TIME = time.time()
# file_ = open('Autocorr_vs_T.xvg', 'w')
filenames = []
for i in range(1, 101):
	filenames.append('proj_%s.xvg' % i)
	# file_.write(str(T) + ' ' + str(R))
	# file_.write('\n')
	# file_.write('\n')
with Pool(8) as p:
	data = p.map(main, filenames)
	for idx, obj in enumerate(data):
		T = obj[0]
		R = obj[1]
		plt.loglog(T, R, color = [0, 0, 1 - idx / 50])
plt.ylim([0.001, 1])
plt.xlabel('Time (ns)')
plt.ylabel('Autocorrelation')
print(time.time() - TIME)
plt.show()