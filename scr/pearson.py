import numpy as np


def main(filename_1, filename_2):
	pearson_value = np.corrcoef(np.loadtxt(filename_1).flatten('F'), np.loadtxt(filename_2).flatten('F'))[0][1]
	print(pearson_value)
	np.savetxt('pearson.dat', np.array([pearson_value]))
	print('Value saved in "pearson.dat"')