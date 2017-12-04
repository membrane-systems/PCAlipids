import numpy as np
import sys

def main(filename_1, filename_2):
	file_1 = open(filename_1, 'r')
	file_2 = open(filename_2, 'r')

	cov_vec_1 = []
	cov_vec_2 = []

	for line in file_1:
		cov_vec_1.extend(line.split())

	for line in file_2:
		cov_vec_2.extend(line.split())	


	for i in range(len(cov_vec_1)):
		cov_vec_1[i] = float(cov_vec_1[i])

	for i in range(len(cov_vec_2)):
		cov_vec_2[i] = float(cov_vec_2[i])

	cov_vec_1 = np.array(cov_vec_1, dtype = np.float64)
	cov_vec_2 = np.array(cov_vec_2, dtype = np.float64)

	print(np.corrcoef(cov_vec_1, cov_vec_2)[0][1])

if __name__ == '__main__':
 	args = sys.argv[1:]
 	if '-f1' in args and '-f2' in args:
 		main(args[args.index('-f1') + 1], args[args.index('-f2') + 1])