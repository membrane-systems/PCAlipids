import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from multiprocessing import Pool
import sys
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
	print(filename + ' - processed')
	# plt.scatter(T[:100],R[:100],s = 5)
	# j = 1
	# while j < len(R[:100]) / 2:
	# 	plt.plot([T[:100][j], T[:100][j*2]], [R[:100][j],R[:100][j * 2]], color = 'red')
	# 	j *= 2
	# plt.yscale('log')
	# plt.xscale('log')
	# plt.show()
	return T, R


def main(filenames, N_lips, timestep, file_out = 'autocorr_relaxtime_vs_PC.xvg'):
	input_data = [(filename, N_lips, timestep) for filename in filenames]
	name = filenames[0][:filenames[0].rfind('_')]
	# print(input_data)
	with Pool(8) as p:
		data = p.starmap(calc, input_data)
	file = open('%s_AUTO_VS_T.xvg' % name, 'w')
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
	plt.savefig('autocorrelation.png')
	plt.clf()

	file = open('%s_' % name + file_out, 'w')
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
	plt.savefig('autocorrelation_relax_time.png')


# if __name__ == '__main__':
# 	args = sys.argv[1:]
# 	if '-p' in args and '-ln' in args and '-dt' in args and '-o' in args and '-pr' not in args:
# 		filenames = [args[i] for i in range(args.index('-p') + 1, args.index('-o'))]
# 		main(filenames, int(args[args.index('-ln') + 1]), float(args[args.index('-dt') + 1]), args[args.index('-o') + 1])
# 	elif '-p' in args and '-ln' in args and '-dt' in args and '-o' not in args and '-pr' not in args:
# 		print('No output file supplied. Data will be written in "autocorr_relaxtime_vs_PC.xvg"')
# 		filenames = [args[i] for i in range(args.index('-p') + 1, args.index('-o'))]
# 		main(filenames, int(args[args.index('-ln') + 1]), float(args[args.index('-dt') + 1]))
# 	elif '-pr' in args and '-ln' in args and '-dt' in args and '-o' in args and '-p' not in args:
# 		files = args[args.index('-pr') + 1]
# 		file_start = files[:files.find('-')]
# 		file_end = files[files.find('-') + 1:]
# 		start = int(''.join(filter(lambda x: x.isdigit(), file_start)))
# 		end = int(''.join(filter(lambda x: x.isdigit(), file_end)))
# 		file_mask = ''.join(filter(lambda x: not x.isdigit(), file_start))
# 		filenames = [file_mask[:file_mask.find('.')] + str(i) + file_mask[file_mask.find('.'):] for i in range(start, end + 1)]
# 		main(filenames, int(args[args.index('-ln') + 1]), float(args[args.index('-dt') + 1]), args[args.index('-o') + 1])
# 	elif '-pr' in args and '-ln' in args and '-dt' in args and '-o' not in args and '-p' not in args: 
# 		files = args[args.index('-pr') + 1]
# 		file_start = files[:files.find('-')]
# 		file_end = files[files.find('-') + 1:]
# 		start = int(''.join(filter(lambda x: x.isdigit(), file_start)))
# 		end = int(''.join(filter(lambda x: x.isdigit(), file_end)))
# 		file_mask = ''.join(filter(lambda x: not x.isdigit(), file_start))
# 		filenames = [file_mask[:file_mask.find('.')] + str(i) + file_mask[file_mask.find('.'):] for i in range(start, end + 1)]
# 		main(filenames, int(args[args.index('-ln') + 1]), float(args[args.index('-dt') + 1]))
# 	elif '-h' not in args:
# 		print('Missing parameters, try -h for flags\n')
# 	else:
# 		print('-p <sequence of projection files>\n -pr <range of files: "proj1.xvg-proj100.xvg">\n-t <topology file> (any file with topology)\n-r <reference traj file> (any file with 1 ref frame). If not supplied, \
# 			the first frame of trajectory will be used for alignment.')

