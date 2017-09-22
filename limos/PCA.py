import mdtraj as md
import numpy as np
import sys
import os
import time


def load_traj(traj_file, top_file):
	traj = md.load(traj_file, top = top_file)
	return traj


def PCA(traj, first_PC, last_PC):
	x_std = traj.xyz.astype(np.float64).reshape(traj.n_frames,traj.n_atoms*3).T
	mean_vec = np.mean(x_std.T,axis=0,dtype=np.float64)
	x_std = x_std-np.array([mean_vec]).T
	cov_mat = np.cov(x_std)
	print("Covariance matrix calculated (%s,%s)" % cov_mat.shape)
	trace = 0.
	for i in range(len(cov_mat)):
		trace += cov_mat[i][i]
	print('Trace of the covariance matrix: %s' % trace)
	eig_vals, eig_vecs = np.linalg.eig(cov_mat)
	if first_PC == None and last_PC == None:
		first_PC = 1
		last_PC = len(eig_vecs)
	elif first_PC == None:
		first_PC = 1
		last_PC = int(last_PC)
	elif last_PC == None:
		first_PC = int(first_PC)
		last_PC = len(eig_vecs)
	else:
		first_PC = int(first_PC)
		last_PC = int(last_PC)
	proj = []
	for i in range(first_PC - 1, last_PC):
		proj.append((x_std).T.dot(eig_vecs[:,i]))
	return proj, eig_vals, eig_vecs, cov_mat, first_PC, last_PC


def main(traj_file, top_file, first_PC = None, last_PC = None):
	PATH = os.getcwd() + '/'
	proj, eig_vals, eig_vecs, cov_mat, first_PC, last_PC = PCA(load_traj(PATH + traj_file, PATH + top_file), first_PC, last_PC)
	with open(PATH + 'eigenval.xvg', 'w') as file:
		file.write('@    title "Eigenvalues of the covariance matrix"\n')
		for i in range(len(eig_vals)):
			file.write(str(i + 1) + ' ' + str(eig_vals[i]) + '\n')
		file.write('&')
	print('Wrote %s eigenvalues in "eigenval.xvg"' % len(eig_vals))
	with open(PATH + 'projection.xvg', 'w') as file:
		file.write('@    Number of projections: %s\n' % len(proj))
		for i in range(len(proj)):
			# print(i)
			file.write('@    title "Projection %s\n' % (first_PC + i))
			file.write(''.join((str(j * 0.1) + '     ' + str(proj[i][j]) + '\n') for j in range(len(proj[i]))))
			# for j in range(len(proj[i])):
			# 	file.write(str(proj[i][j]) + '\n')
			file.write('&\n')
	print('Wrote %s projections in "projection.xvg"' % len(proj))
	with open(PATH + 'covar.dat', 'w') as file:
		flag = 0
		for i in range(len(cov_mat)):
			for j in range(len(cov_mat)):
				file.write(str(cov_mat[i][j]))
				flag += 1
				if flag == 3:
					flag = 0
					file.write('\n')
	print('Wrote covariance matrix in "covar.dat"')
	return 'projection.xvg', first_PC, last_PC


# if __name__ == '__main__':
# 	# start = time.time()
# 	args = sys.argv[1:]
# 	if '-f' in args and '-t' in args and '-first' in args and '-last' in args:
# 		main(args[args.index('-f') + 1], args[args.index('-t') + 1], args[args.index('-first') + 1], args[args.index('-last') + 1])
# 	elif '-f' in args and '-t' in args and '-first' not in args and '-last' not in args:
# 		main(args[args.index('-f') + 1], args[args.index('-t') + 1], None, None)
# 	elif '-f' in args and '-t' in args and '-first' in args and '-last' not in args:
# 		main(args[args.index('-f') + 1], args[args.index('-t') + 1], args[args.index('-first') + 1], None)
# 	elif '-f' in args and '-t' in args and '-first' not in args and '-last' in args:
# 		main(args[args.index('-f') + 1], args[args.index('-t') + 1], None, args[args.index('-last') + 1])
# 	elif '-h' not in args:
# 		print('Missing parameters, try -h for flags\n')
# 	else:
# 		print('-f <trajectory file> (file format *.xtc)\n-t <topology file> (any file with topology)\n -first <first PC> -last <last PC>.')