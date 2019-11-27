import numpy as np
import math


def get_data_from_file(file_name):
	data = []  
	file_ = open(file_name, 'r')
	lines = file_.readlines()
	file_.close()
	n = 0
	for i in range(2, len(lines)):
		n += 1
		if lines[i].find('*') == -1 and lines[i] != '\n':
			data.append(float(lines[i].split()[1]))
		else:
			break
	return data[:n]


def get_data_from_file_1(file_name):
	data = []  
	file_ = open(file_name, 'r')
	lines = file_.readlines()[1:-1]
	file_.close()
	n = 0
	for line in lines:
		n += 1
		data.append(float(line.split()[1]))
	return data[:n]



def f(data_1, data_2, eigvals,outF):
	beta = 0.
	
	m = min(len(data_1),len(data_2),len(eigvals))
	if (m < len(eigvals)):
		print("It is recommended to calculate the relative timescales for all projections")
	for eigval in eigvals[:m]:
		beta += eigval
	beta = 1. / beta
	R = 0.
	for i in range(m):
		R += beta * eigvals[i] * math.log(data_1[i]/data_2[i])
	file = open(outF, 'w')
	file.write(str(math.exp(R)) + ' ' + str(1 / math.exp(R)))
	file.close()
	return math.exp(R)


def main(file_eigvals, file1, file2,outF):
	if not file_eigvals or not file1 or not file2:
		print("Eigenvalues and characteristic timescales for 2 trajectories have to be provided.\n\
Run pcalipids.py reltime -h for help")
		return 0
	data_1 = get_data_from_file(file1)
	data_2 = get_data_from_file(file2)
	eigvals = get_data_from_file_1(file_eigvals)
	print('Comparison value : ' + str(f(data_1, data_2, eigvals,outF)))
	print('Values saved in '+outF)
