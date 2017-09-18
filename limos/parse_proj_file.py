import numpy as np
from scipy.stats.kde import gaussian_kde
import matplotlib.pyplot as plt
from multiprocessing import Pool
import sys
import os


def get_data_from_file(file_name, proj_num):
	file = open(file_name, 'r')
	line = file.readline()
	i = 0
	projs = []
	while i < proj_num:
		proj = []
		while line.find('&') != 0:
			if line.find('@') != 0:
				proj.append(float(line.split()[1]))
			line = file.readline()
		proj.append(i + 1)
		i += 1
		line = file.readline()
		projs.append(proj)
	file.close()
	return projs


def plot_dist(data, PATH):
	N = data[len(data) - 1]
	data = np.array(data[:len(data) - 1])
	KDEpf = gaussian_kde(data)
	x = np.linspace(np.amin(data) - 0.5, np.amax(data) + 0.5, 200)
	plt.ylabel('Probability density (a.u.)')
	plt.xlabel('PC projection value (A')
	plt.title('Distributions of the PC%s projection' % N)
	plt.plot(x, KDEpf(x), 'r', color = 'blue')
	plt.savefig(PATH + 'PC%s_dist.png' % N)
	plt.clf()


def main(file_name, proj_num):
	PATH = os.getcwd() + '/'
	projs = get_data_from_file(PATH + file_name, int(proj_num))
	# for data in projs:
	# 	plot_dist(data, PATH)
	for i in range(len(projs)):
		projs[i] = (projs[i], PATH)
	with Pool(2) as p:
		p.starmap(plot_dist, projs)


if __name__ == '__main__':
	args = sys.argv[1:]
	if '-p' in args and '-n' in args:
		main(args[args.index('-p') + 1], args[args.index('-n') + 1])
	elif '-p' in args and '-n' not in args:
		main(args[args.index('-p') + 1], 3)
	elif '-h' not in args:
		print('Missing parameters, try -h for flags\n')
	else:
		print('-p <projection file> (file format *.xvg)\n-h <number of first projection> (int format). If not supplied, the firts 3 projections will be analyzed.')
	