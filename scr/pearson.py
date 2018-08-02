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
	file = open('pearson.dat', 'w')
	file.write(str(np.corrcoef(cov_vec_1, cov_vec_2)[0][1]))
	file.close()
	print('Value saved in "pearson.dat"')
	
# if __name__ == '__main__':
#  	args = sys.argv[1:]
#  	if '-cov1' in args and '-cov2' in args:
#  		main(args[args.index('-cov1') + 1], args[args.index('-cov2') + 1])
#  	elif '-h' not in args and '-help' not in args:
#  		print('Missing parameters, try -h for flags\n')
#  	else:
#  		print('-cov1, -cov2 - 2 files with covariance matrices')