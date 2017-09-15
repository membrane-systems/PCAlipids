import numpy
import matplotlib.pyplot as plt
from scipy.stats.kde import gaussian_kde
import time 
import sys
import numpy as np


def get_data_from_file(file_name):
	file = open(file_name, 'r')
	lines = file.readlines()
	file.close()
	lines = lines[:len(lines) - 1]
	while '' in lines:
		lines.remove()
	data = []
	for line in lines:
		data.append(float(line))
	return np.array(data)


def PDFs(data):
	KDEpdf = gaussian_kde(data)
	max_ = np.amax(data)
	min_ = np.amin(data)
	print(max_)
	y = np.linspace(min_ - 0.5, max_ + 0.5, 3000)
	plt.plot(y, KDEpdf(y), 'r', label = "PC1 dist", color = "blue")
	plt.show()


def main(file_name):
	PDFs(get_data_from_file(file_name))


if __name__ == "__main__":
	args = sys.argv[1:]
	if '-p' in args:
		main(args[args.index('-p') + 1])
	elif '-h' not in args:
		print('Missing parameters, try -h for flags\n')
	else:
		print('-p <projection file> (file format *.xvg)')