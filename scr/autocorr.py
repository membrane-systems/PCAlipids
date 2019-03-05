import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from multiprocessing import Pool
import sys, os
import math


def estimated_autocorrelation(x, variance, mean):
	n = len(x)
	# variance = x.var()
	# print(variance)
	x = x-mean
	r = signal.fftconvolve(x, x[::-1], mode="full")[-n:]
	# assert np.allclose(r, np.array([(x[:n-k]*x[-(n-k):]).sum() for k in range(n)]))
	result = r/(variance*(np.arange(n, 0, -1)))
	return result


def get_nearest_value(iterable, value):
	# idx = (np.abs(iterable - value)).argmin()
	# B = idx
	# if B == 0:
	# 	A = 1
	# elif B == len(iterable) - 1:
	# 	A = B - 1
	# elif (iterable[idx] - value) * (iterable[idx + 1] - value) > 0:
	# 	A = idx - 1
	# else:
	# 	A = idx + 1
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
	data = []
	for line in file:
	    if line.find('@') == -1 and line.find('&') == -1:
	        data.append(np.float(line.split()[0]))
	    if line.find('&') != -1:
	        break
	file.close()
	data = np.array(data, dtype = np.float64)
	variance  = data.var()
	mean = data.mean()
	data_1 = data[:len(data) // N_lips]
	R = estimated_autocorrelation(data_1, variance, mean)
	for i in range(1, N_lips):
	    data_1 = data[len(data) * i // N_lips:len(data) * (i + 1) // N_lips]
	    R += estimated_autocorrelation(data_1, variance, mean)
	R /= N_lips
	# print(len(R))
	T = []
	for i in range(0, len(R)):
	    T.append(timestep * i)
	print(filename.split('/')[-1] + ' - processed')
	# plt.scatter(T[:100],R[:100],s = 5)
	# j = 1
	# while j < len(R[:100]) / 2:
	# 	plt.plot([T[:100][j], T[:100][j*2]], [R[:100][j],R[:100][j * 2]], color = 'red')
	# 	j *= 2
	# plt.yscale('log')
	# plt.xscale('log')
	# plt.show()
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

	with Pool(8) as p:
		data = p.starmap(calc, input_data)

	file = open(name+'_values_vs_t'+rez, 'w')
	for j in range(len(data[0][0])):
		for i in range(len(data)):
			if i == 0:
				file.write(str(data[i][0][j]) + ' ' + str(str(data[i][1][j])))
			else:
				file.write(' ' + str(str(data[i][1][j])))
		file.write('\n')
	file.close()
	POINTS = []
	j = 1
	while j < len(data[0][0]):
		POINTS.append(j)
		j = int(1.5 * j) + 1
	# print(POINTS)
	for idx, obj in enumerate(data[:10]):
		T = obj[0]
		R = obj[1]
		# T_pic = [T[i] for i in POINTS]
		# R_pic = [R[i] for i in POINTS]
		plt.loglog(T, R, color = [0, 0, 1 - idx / 10])
	plt.ylim([0.001, 1])
	plt.xlabel('Time (ns)')
	plt.ylabel('Autocorrelation')
	plt.savefig(name+'_values_vs_t'+'.png')
	plt.clf()

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
		A,B = get_nearest_value(R_log, -2)
		a = (T_log[B] - T_log[A]) / (R_log[B] - R_log[A])
		# print(T_log[B], T_log[A], R_log[B], R_log[A])
		# print(A)
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
		# print(math.e ** t_relax)
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
