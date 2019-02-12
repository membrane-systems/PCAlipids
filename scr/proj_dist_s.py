import numpy as np
from scipy.stats.kde import gaussian_kde
import matplotlib.pyplot as plt
from multiprocessing import Pool
import sys
import os


def get_data_from_file(file_name):
	file = open(file_name, 'r')
	proj = []
	line = file.readline()
	while line.find('&') == -1:
		if line.find('@') == -1:
			proj.append(float(line.split()[0]))
		line = file.readline()
	file.close()
	return proj


def plot_dist(data, PATH, PC):
	N = PC
	data = np.array(data)
	KDEpf = gaussian_kde(data)
	x = np.linspace(np.amin(data), np.amax(data), 100)
	plt.ylabel('Probability density (a.u.)')
	plt.xlabel('PC projection value (A)')
	plt.title('Distribution of the PC%s projection' % N)
	plt.plot(x, KDEpf(x), 'r', color = 'blue')
	plt.savefig(PATH + 'PC%s_dist.png' % N)
	plt.clf()


def main(file_name, PC):
	PATH = os.getcwd() + '/'
	proj = get_data_from_file(PATH + file_name)
	proj = [(proj,PATH,PC)]
	with Pool(2) as p:
		p.starmap(plot_dist, proj)
