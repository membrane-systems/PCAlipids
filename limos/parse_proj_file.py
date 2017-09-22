import numpy as np
from scipy.stats.kde import gaussian_kde
import matplotlib.pyplot as plt
from multiprocessing import Pool
import sys
import os


def get_data_from_file(file_name, first_PC, last_PC):
	file = open(file_name, 'r')
	line = file.readline()
	i = 1
	projs = []
	while i < first_PC:
		line = file.readline()
		if line.find('&') != -1:
			i += 1
	line = file.readline()
	while i <= last_PC:
		proj = []
		while line.find('&') == -1:
			if line.find('@') == -1:
				proj.append(float(line.split()[1]))
			line = file.readline()
		proj.append(i)
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
	plt.xlabel('PC projection value (A)')
	plt.title('Distribution of the PC%s projection' % N)
	plt.plot(x, KDEpf(x), 'r', color = 'blue')
	plt.savefig(PATH + 'PC%s_dist.png' % N)
	plt.clf()


def main(file_name, first_PC, last_PC):
	PATH = os.getcwd() + '/'
	projs = get_data_from_file(PATH + file_name, int(first_PC), int(last_PC))
	for i in range(len(projs)):
		projs[i] = (projs[i], PATH)
	with Pool(2) as p:
		p.starmap(plot_dist, projs)


# if __name__ == '__main__':
# 	args = sys.argv[1:]
# 	if '-p' in args and '-first' in args and '-last' in args:
# 		main(args[args.index('-p') + 1], args[args.index('-first') + 1], args[args.index('-last') + 1])
# 	elif '-p' in args and '-first' not in args and '-last' not in args:
# 		main(args[args.index('-p') + 1], 1, 3)
# 	elif '-h' not in args:
# 		print('Missing parameters, try -h for flags\n')
# 	else:
# 		print('-p <projection file> (file format *.xvg)\n -fisrt <first projection> -last <last projection> (int format). \nIf not supplied, the first 3 projections will be analyzed.')