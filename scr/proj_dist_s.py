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


def main(file_name, range_fnames):
	PATH = os.getcwd() + '/'
	# if both file or range of files are not defined 
	if not file_name and not range_fnames:
		print("The projections files have to be provided.\n\
		Use pcalipids.py projdist -h for help")
	# if range is defined
	elif range_fnames:
		# define file name for the first projection file in list
		firstFile = range_fnames[:range_fnames.find('-')]
		# define file name for the last projection file in list
		lastFile  = range_fnames[range_fnames.find('-')+1:]
		# find first PC
		first = int(firstFile[firstFile.rfind('_')+1:firstFile.rfind('.')])
		# find last PC
		last = int(lastFile[lastFile.rfind('_')+1:lastFile.rfind('.')])
		# define filemask
		mask = firstFile[:firstFile.rfind('_')]
		# define filerez
		rez = firstFile[firstFile.rfind('.'):]
		for i in range(first,last+1):
			proj = get_data_from_file(PATH + mask + "_" + str(i) + rez)
			proj = [(proj,PATH,i)]
			with Pool(2) as p:
				p.starmap(plot_dist, proj)
	# if file is defined
	else:
		proj = get_data_from_file(PATH + file_name)
		PC = int(file_name[file_name.rfind('_')+1:file_name.rfind('.')])
		proj = [(proj,PATH,PC)]
		with Pool(2) as p:
			p.starmap(plot_dist, proj)
