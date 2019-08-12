import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from multiprocessing import Pool
import sys, os
import math


def estimated_autocorrelation(x, variance, mean):
	# x - data
	n = len(x)
	# Center the data
	x = x-mean
	# Convolve the data
	r = signal.fftconvolve(x, x[::-1], mode="full")[-n:]
	# Weight the correlation to get the result in range -1 to 1
	result = r/(variance*(np.arange(n, 0, -1)))
	return result

# Interpolation
def get_nearest_value(iterable, value):
	for idx, x in enumerate(iterable):
		if x < value:
			break

	A = idx
	if idx == 0:
		for idx, x in enumerate(iterable):
			if x < iterable[A]:
				break
		B = idx
	else:
		B = idx - 1
	# print(A,B)
	return A, B


def calc(filename, N_lips, timestep):
	file = open(filename, 'r')

	# Colecting the data
	data = []
	for line in file:
	    if line.find('@') == -1 and line.find('&') == -1:
	        data.append(np.float(line.split()[0]))
	    if line.find('&') != -1:
	        break
	file.close()
	data = np.array(data, dtype = np.float64)
	
	# Calculate data parameters
	variance  = data.var()
	mean = data.mean()
	
	# Get correlation for the first lipid
	data_1 = data[:len(data) // N_lips]
	R = estimated_autocorrelation(data_1, variance, mean)
	
	# Get correlation for other lipids
	for i in range(1, N_lips):
	    data_1 = data[len(data) * i // N_lips:len(data) * (i + 1) // N_lips]
	    R += estimated_autocorrelation(data_1, variance, mean)
	
	# Calculate averaged correlation
	R /= N_lips
	
	# Get corresponding times
	T = []
	for i in range(0, len(R)):
	    T.append(timestep * i)
	print(filename.split('/')[-1] + ' - processed')
	return T, R


def main(file_name, range_fnames, N_lips, timestep, file_out ):
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
		input_data = [(PATH + mask + "_" + str(i) + rez, \
			N_lips, timestep) \
			for i in range(first,last+1)]
	else:
		input_data = [(file_name, N_lips, timestep)]

	name = file_out[:file_out.rfind('.')]
	rez  = file_out[file_out.rfind('.'):]

	# Calculate correlation
	with Pool(8) as p:
		data = p.starmap(calc, input_data)

	# Write correlation to file
	file = open(name+'_values_vs_t'+rez, 'w')
	for j in range(len(data[0][0])):
		for i in range(len(data)):
			if i == 0:
				file.write(str(data[i][0][j]) + ' ' + str(str(data[i][1][j])))
			else:
				file.write(' ' + str(str(data[i][1][j])))
		file.write('\n')
	file.close()
	
	# Interpolation of data
	# Create time array
	POINTS = []
	j = 1
	while j < len(data[0][0]):
		POINTS.append(j)
		j = int(1.5 * j) + 1

	# Plotting autocorrelation
	for idx, obj in enumerate(data[:10]):
		T = obj[0]
		R = obj[1]
		plt.loglog(T, R, color = [0, 0, 1 - idx / 10])

	plt.ylim([0.001, 1])
	plt.xlabel('Time (ns)')
	plt.ylabel('Autocorrelation')
	plt.savefig(name+'_values_vs_t'+'.png')
	plt.clf()

	# Save autocorrelation decay times to file
	# Write data for e^2 decay
	file = open(name + '_relax_times_vs_pc' + rez, 'w')
	file.write('### Autocorr ###\n')
	file.write('E**2\n')

	handles = []

	PC = []
	T_relax = []

	for i, value in enumerate(data):
		T = value[0]
		R = value[1]
		T_pic = np.array([T[i] for i in POINTS if R[i] > 0.])
		R_pic = np.array([R[i] for i in POINTS if R[i] > 0.])
		R_log = np.log(R_pic[:])
		T_log = np.log(T_pic[:])
		# Interpolate data
		A,B = get_nearest_value(R_log, -2)
		a = (T_log[B] - T_log[A]) / (R_log[B] - R_log[A])
		b = T_log[A] - a * R_log[A]
		t_relax = a * -2 + b
		PC.append(i + 1)
		T_relax.append(math.e ** t_relax)

	p = plt.loglog(PC[:10], T_relax[:10], label = r'$\tau_2 = 1/e^2$', color = 'blue', linestyle = '--')
	handle, = p 
	handles.append(handle)

	for i in range(len(PC)):
		file.write(str(PC[i]) + ' ' + str(T_relax[i]))
		file.write('\n')
	file.write('\n')

	# Write data for e decay
	PC = []
	T_relax = []

	for i, value in enumerate(data):
		T = value[0]
		R = value[1]
		T_pic = np.array([T[i] for i in POINTS if R[i] > 0.])
		R_pic = np.array([R[i] for i in POINTS if R[i] > 0.])
		R_log = np.log(R_pic[:])
		T_log = np.log(T_pic[:])
		A,B = get_nearest_value(R_log,-1)
		a = (T_log[B] - T_log[A]) / (R_log[B] - R_log[A])
		b = T_log[A] - a * R_log[A]
		t_relax = a * -1 + b
		PC.append(i + 1)
		T_relax.append(math.e ** t_relax)


	file.write('E**1\n')
	for i in range(len(PC)):
		file.write(str(PC[i]) + ' ' + str(T_relax[i]))
		file.write('\n')
	file.close()


	p = plt.loglog(PC[:10], T_relax[:10], label = r'$\tau_1 = 1/e$', color = 'blue', linestyle = '-')
	handle, = p
	handles.append(handle)
	plt.ylim([0.01, 10 ** 3])
	plt.ylabel('Relaxation time (ns)')
	plt.xlabel('Component')
	plt.legend(handles = handles)
	plt.savefig(name + '_relax_times_vs_pc.png')
