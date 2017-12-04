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


def main(traj_file, top_file, val_file, vec_file, cov_file, proj_file, first_PC = None, last_PC = None):
	PATH = os.getcwd() + '/'
	proj, eig_vals, eig_vecs, cov_mat, first_PC, last_PC = PCA(load_traj(PATH + traj_file, PATH + top_file), first_PC, last_PC)
	
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

	if proj_file == None:
		file_out = 'projection.xvg'
	else:
		file_out = proj_file
	with open(PATH + file_out, 'w') as file:
		file.write('@    Number of projections: %s\n' % len(proj))
		for i in range(len(proj)):
			# print(i)
			file.write('@    title "Projection %s\n' % (first_PC + i))
			file.write(''.join((str(j * 0.1) + '     ' + str(proj[i][j]) + '\n') for j in range(len(proj[i]))))
			# for j in range(len(proj[i])):
			# 	file.write(str(proj[i][j]) + '\n')
			file.write('&\n')
	print('Wrote %s projections in "%s"' % (len(proj),file_out))

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
	with open(PATH + file_out, 'w') as file:
		for i in range(len(eig_vecs)):
			for j in range(len(eig_vecs[:,i])):
				file.write(str(eig_vecs[:,i][j]) + ' ')
			file.write('\n')
	print('Wrote eigenvectors in "%s"' % file_out)
	
	return 


if __name__ == '__main__':
	args = sys.argv[1:]
	if '-oeval' in args:
		val_file = args[args.index('-oeval') + 1]
	else:
		val_file = None
	if '-oevec' in args:
		vec_file = args[args.index('-oevec') + 1]
	else:
		vec_file = None
	if '-ocov' in args:
		cov_file = args[args.index('-ocov') + 1]
	else:
		cov_file = None
	if '-op' in args:
		proj_file = args[args.index('-op') + 1]
	else:
		proj_file = None
	if '-f' in args and '-t' in args and '-first' in args and '-last' in args:
		main(args[args.index('-f') + 1], args[args.index('-t') + 1], first_PC = args[args.index('-first') + 1], last_PC = args[args.index('-last') + 1], val_file = val_file, vec_file = vec_file, cov_file = cov_file, proj_file = proj_file)
	elif '-f' in args and '-t' in args and '-first' not in args and '-last' not in args:
		main(args[args.index('-f') + 1], args[args.index('-t') + 1], val_file = val_file, vec_file = vec_file, cov_file = cov_file, proj_file = proj_file, first_PC = None, last_PC = None)
	elif '-f' in args and '-t' in args and '-first' in args and '-last' not in args:
		main(args[args.index('-f') + 1], args[args.index('-t') + 1], first_PC = args[args.index('-first') + 1], val_file = val_file, vec_file = vec_file, cov_file = cov_file, proj_file = proj_file, last_PC = None)
	elif '-f' in args and '-t' in args and '-first' not in args and '-last' in args:
		main(args[args.index('-f') + 1], args[args.index('-t') + 1], last_PC = args[args.index('-last') + 1], val_file = val_file, vec_file = vec_file, cov_file = cov_file, proj_file = proj_file, first_PC =None)
	elif '-h' not in args:
		print('Missing parameters, try -h for flags\n')
	else:
		print('-f <trajectory file> (file format *.xtc)\n-t <topology file> (any file with topology)\n -first <first PC> -last <last PC> \n -oeval <output file with eigenvalues>\n\
			-oevec <output file with eigenvectors>\n -ocov <output file with covariance matrix>\n -op <output file with projections>.')
