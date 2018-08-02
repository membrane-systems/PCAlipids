import numpy as np
import math


def get_data_from_file(file_name):
  data = []  
  file_ = open(file_name, 'r')
  lines = file_.readlines()
  file_.close()
  for i in range(2, len(lines)):
    if lines[i].find('*') == -1 and lines[i] != '\n':
      data.append(float(lines[i].split()[1]))
    else:
      break
  return data[:162]


def get_data_from_file_1(file_name):
  data = []  
  file_ = open(file_name, 'r')
  lines = file_.readlines()[1:-1]
  file_.close()
  for line in lines:
    data.append(float(line.split()[1]))
  return data[:162]



def f(data_1, data_2, eigvals):
  beta = 0.
  for eigval in eigvals:
    beta += eigval
  beta = 1. / beta
  R = 0.
  for i in range(162):
    R += beta * eigvals[i] * math.log(data_1[i]/data_2[i])
  file = open('timescales_comparison.dat', 'w')
  file.write(str(math.exp(R)) + ' ' + str(1 / math.exp(R)))
  file.close()
  return math.exp(R)


def main(file_eigvals, file1, file2):
  data_1 = get_data_from_file(file1)
  data_2 = get_data_from_file(file2)
  eigvals = get_data_from_file_1(file_eigvals)
  print('Comparison value : ' + str(f(data_1, data_2, eigvals)))
  print('Values saved in "timescales_comparison.dat"')