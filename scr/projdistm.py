import numpy
import matplotlib.pyplot as plt
from scipy.stats.kde import gaussian_kde
import time 
import sys
import numpy as np
import matplotlib.patches as mpatches
import math


def get_data_from_file(file_name):
	file = open(file_name, 'r')
	# lines = file.readlines()
	# file.close()
	# flag = False
	# lines = lines[:len(lines) - 1]
	# while '' in lines:
	# 	lines.remove()
	# data = []
	# for line in lines:
	# 	if line.find('@') == -1 and line.find('&') == -1:
	# 		print(line.split())

	# 		data.append(float(line.split()[1]))
	# 	if line.find('&') != -1:
	# 		break

	# return np.array(data)

	data = []
	for line in file:
		data.append(np.float(line))
	file.close()
	return np.array(data)


def PDFs(y, data, label, col, style = '-'):
	KDEpdf = gaussian_kde(data)
	max_ = np.amax(data)
	min_ = np.amin(data)
	K = KDEpdf(y)
	p= plt.plot( - y, KDEpdf(y), 'r', label = label, color = col) 
	return p

def main(file_seq):
	y = np.linspace(-5, 5, 51)
	handles = []
	for file in file_seq:
		line1, = PDFs(y, get_data_from_file(file),file[:file.find('.')], (np.random.random_sample(),np.random.random_sample(), np.random.random_sample()), style = '-')
		handles.append(line1)
		
	plt.xlabel('PC projection value (A)')
	plt.ylabel('Probability density (a.u.)')
	plt.ylim([0,0.63])
	plt.xlim([-4,4.1])
	plt.legend(handles = handles, ncol = 2)
	plt.show()


# if __name__ == "__main__":
# 	args = sys.argv[1:]
# 	main()