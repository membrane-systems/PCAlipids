import numpy as np
import sys

def main(filename_1, filename_2,ofile):
	if not filename_1 or not filename_2:
		print("Two covariance matrices have to be provided.\n\
Run pcalipids.py pearson -h for help")
		return 0
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
	file = open(ofile, 'w')
	file.write(str(np.corrcoef(cov_vec_1, cov_vec_2)[0][1]))
	file.close()
	print('Value saved in '+ofile)
	
