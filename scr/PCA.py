import mdtraj as md
import numpy as np
import sys
import os
import time

# Rudiment function - not using now
def PCA(traj):
	x_std = traj.xyz.astype(np.float64).reshape(traj.n_frames,traj.n_atoms*3).T
	mean_vec = np.mean(x_std.T,axis=0,dtype=np.float64)
	x_std = x_std-np.array([mean_vec]).T
	cov_mat = np.cov(x_std)
	print("Covariance matrix calculated (%s,%s)" % cov_mat.shape)
	trace = np.matrix.trace(cov_mat)
	print('Trace of the covariance matrix: %s' % trace)
	eig_vals, eig_vecs = np.linalg.eigh(cov_mat)
	return eig_vals[::-1], eig_vecs, cov_mat


def PCA_mem(traj, top):
	# get coordinates of average structure
	mean_str = md.load(top)
	mean_vec = mean_str.xyz.astype(np.float64).\
	reshape(1, mean_str.n_atoms * 3)
	
	firstLoad = True # check if we load frames for the first times
	N = 0 # number of loaded frames

	for frame in md.iterload(traj, top = top, chunk = 100000):
		N += frame.n_frames
		# when we start loading trajectory we need to initialize
		# mean vectors and X.(X^T) matrices
		if firstLoad:
			X_1 = np.array([0.0] * frame.n_atoms * 3,dtype = np.float64)
			X_X = np.array([[0.0] * frame.n_atoms * 3 for i in \
				range(frame.n_atoms * 3)], dtype = np.float64)
			firstLoad = False

		# X - coordinates of atoms subjacted average structure
		# X_1 - the sum of all coordinates (to calculate mean)
		# X_X - product of X and X^T
		X = frame.xyz.astype(np.float64).reshape\
		(frame.n_frames, frame.n_atoms * 3) - mean_vec
		X_1 += X.sum(axis=0)
		X_X += np.tensordot(X, X, axes=(0,0))

	# covariance matrix calculation
	cov_mat = np.empty((len(X_1), len(X_1)), dtype = np.float64)
	cov_mat = (X_X - np.dot(X_1.reshape(len(X_1),1), \
		(X_1.reshape(len(X_1),1)).T) / N) / (N - 1)
	print("Covariance matrix calculated (%s,%s)" % cov_mat.shape)
	
	trace = np.matrix.trace(cov_mat)
	print('Trace of the covariance matrix: %s' % trace)
	
	# calculation of eigenvalues and eigenvectors
	eig_vals, eig_vecs = np.linalg.eigh(cov_mat)
	return eig_vals[::-1], eig_vecs, cov_mat


def main(traj_file, top_file, val_file, vec_file, cov_file):
	if not traj_file or not top_file:
		print("Trajectory and average structure have to be provided.\n\
Run pcalipids.py covar -h for help")
		return 0

	PATH = os.getcwd() + '/'
	eig_vals, eig_vecs, cov_mat = PCA_mem(PATH + traj_file, PATH + top_file)
	
	# save eigenvalues	
	np.savetxt(PATH + val_file,\
		np.array([np.arange(1,len(eig_vals)+1),eig_vals]).T,\
		header='@    title "Eigenvalues of the covariance matrix"',\
		footer='&', comments='',fmt=['%5u','%20.17g'])
	print('Wrote %s eigenvalues in "%s"' % (len(eig_vals), val_file))

	# save covariance matrix
	np.savetxt(PATH + cov_file,\
		cov_mat.reshape(cov_mat.shape[0]*cov_mat.shape[0]//3,3), \
		fmt=['%20.17g','%20.17g','%20.17g'])
	print('Wrote covariance matrix in "%s"' % cov_file)

	# save eigenvectors
	np.savetxt(PATH + vec_file,\
		np.flip(eig_vecs,axis=1).T,fmt='%20.17g')
	print('Wrote eigenvectors in "%s"' % vec_file)
	return 
