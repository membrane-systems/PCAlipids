import numpy as np
import sys
import math as m
import scipy as sp

def main(filename_1, filename_2,ofile):
	if not filename_1 or not filename_2:
		print("Two covariance matrices have to be provided.\n\
Run pcalipids.py pearson -h for help")
		return 0
	
	cov1 = np.loadtxt(filename_1)
	cov2 = np.loadtxt(filename_2)

	ndim = int(m.sqrt(cov1.shape[0]*cov1.shape[1]))

	cov1=cov1.reshape((ndim,ndim))
	cov2=cov2.reshape((ndim,ndim))

	sq1 = np.array(sp.linalg.sqrtm(cov1))
	sq2 = np.array(sp.linalg.sqrtm(cov2))

	d = np.trace(np.linalg.matrix_power(sq1-sq2,2))
	overlap = 1. - d/(m.sqrt(np.trace(cov1)+np.trace(cov2)))

	print(overlap)
	np.savetxt(ofile,np.array([overlap]))

	print('Value saved in '+ofile)
	
