import mdtraj as md
import numpy as np
import sys
import os
import time


def load_traj(traj_file, top_file):
	traj = md.load(traj_file, top = top_file)
	return traj


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
	print("The PCA is performed while saving memory, it may take some time.")
	flag = False
	N = 0
	mean_str = md.load(top)
	mean_vec = mean_str.xyz.astype(np.float64).reshape(1, mean_str.n_atoms * 3)
	for frame in md.iterload(traj, top = top, chunk = 100000):
		N += frame.n_frames
		if not flag:
			X_1 = np.array([0.0] * frame.n_atoms * 3,dtype = np.float64)
			X_X = np.array([[0.0] * frame.n_atoms * 3 for i in range(frame.n_atoms * 3)],dtype = np.float64)
			flag = True
		X = frame.xyz.astype(np.float64).reshape(frame.n_frames, frame.n_atoms * 3) - mean_vec
		X_1 += X.sum(axis=0)
		X_X += np.tensordot(X, X, axes=(0,0))
	cov_mat = np.empty((len(X_1), len(X_1)), dtype = np.float64)
	cov_mat = (X_X - np.dot(X_1.reshape(len(X_1),1), (X_1.reshape(len(X_1),1)).T) / N) / (N - 1)
	print("Covariance matrix calculated (%s,%s)" % cov_mat.shape)
	trace = np.matrix.trace(cov_mat)
	print('Trace of the covariance matrix: %s' % trace)
	eig_vals, eig_vecs = np.linalg.eigh(cov_mat)
	return eig_vals[::-1], eig_vecs, cov_mat


def main(traj_file, top_file, val_file, vec_file, cov_file, invert, memory_flag):
	PATH = os.getcwd() + '/'
	if memory_flag == False:
		eig_vals, eig_vecs, cov_mat = PCA(load_traj(PATH + traj_file, PATH + top_file))
	else:
		eig_vals, eig_vecs, cov_mat = PCA_mem(PATH + traj_file, PATH + top_file)
	if val_file == None:
		file_out = 'eigenval.xvg'
	else:
		file_out = val_file
	with open(PATH + file_out, 'w') as file:
		file.write('@    title "Eigenvalues of the covariance matrix"\n')
		for i in range(len(eig_vals)):
			file.write(str(i + 1) + ' ' + str(eig_vals[i]) + '\n')
		file.write('&')
	print('Wrote %s eigenvalues in "%s"' % (len(eig_vals), file_out))

	if cov_file == None:
		file_out = 'covar.dat'
	else:
		file_out = cov_file
	with open(PATH + file_out, 'w') as file:
		flag = 0
		for i in range(len(cov_mat)):
			for j in range(len(cov_mat)):
				file.write(str(cov_mat[i][j]) + ' ')
				flag += 1
				if flag == 3:
					flag = 0
					file.write('\n')
	print('Wrote covariance matrix in "%s"' % file_out)

	if vec_file == None:
		file_out = 'eigenvec.xvg'
	else:
		file_out = vec_file
	if invert:
		invert = -1.
	else:
		invert = 1.
	with open(PATH + file_out, 'w') as file:
		for i in range(len(eig_vecs)):
			for j in range(len(eig_vecs[:,len(eig_vecs) - 1 - i])):
				file.write(str(invert * eig_vecs[:,len(eig_vecs) - 1 - i][j]) + ' ')
			file.write('\n')
	print('Wrote eigenvectors in "%s"' % file_out)
	return 


# if __name__ == '__main__':
# 	args = sys.argv[1:]
# 	if '-oeval' in args:
# 		val_file = args[args.index('-oeval') + 1]
# 	else:
# 		val_file = None
# 	if '-oevec' in args:
# 		vec_file = args[args.index('-oevec') + 1]
# 	else:
# 		vec_file = None
# 	if '-ocov' in args:
# 		cov_file = args[args.index('-ocov') + 1]
# 	else:
# 		cov_file = None
# 	if '-f' in args and '-t' in args:
# 		main(args[args.index('-f') + 1], args[args.index('-t') + 1], val_file = val_file, vec_file = vec_file, cov_file = cov_file)
# 	elif '-h' not in args:
# 		print('Missing parameters, try -h for flags\n')
# 	else:
# 		print('-f <trajectory file> (file format *.xtc)\n-t <topology file> (any file with topology)\n -first <first PC> -last <last PC> \n -oeval <output file with eigenvalues>\n\
# 			-oevec <output file with eigenvectors>\n -ocov <output file with covariance matrix>\n -op <output file with projections>.')

