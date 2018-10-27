import mdtraj as md
import numpy as np
import sys
import os


def load_traj(traj_file, top_file):
	traj = md.load(traj_file, top = top_file)
	return traj


def load_evecs(evecs):
	eigenvecs = []
	file = open(evecs, 'r')
	for line in file:
		a = line.split()
		for i in range(len(a)):
			a[i] = np.float64(a[i])
		a = np.array(a, dtype = np.float64)
		eigenvecs.append(a)
	file.close()
	eigenvecs = np.array(eigenvecs, dtype = np.float64)
	return eigenvecs


def get_proj(traj, first_PC, last_PC, aver, evecs):
	ref_aver_str = md.load(aver)
	traj = traj.superpose(ref_aver_str)
	mean_vec = ref_aver_str.xyz.astype(np.float64)
	x_std = traj.xyz.astype(np.float64).reshape(traj.n_frames,traj.n_atoms*3).T
	x_std = x_std - np.array([mean_vec.reshape(traj.n_atoms*3, )]).T
	eigenvecs = load_evecs(evecs)
	if first_PC == None and last_PC == None:
		first_PC = 1
		last_PC = len(eigenvecs)
	elif first_PC == None:
		first_PC = 1
		last_PC = int(last_PC)
	elif last_PC == None:
		first_PC = int(first_PC)
		last_PC = len(eigenvecs)
	else:
		first_PC = int(first_PC)
		last_PC = int(last_PC)
	proj = []
	for i in range(first_PC - 1, last_PC):
		proj.append((x_std).T.dot(eigenvecs[i]))
	return proj, first_PC, last_PC


def get_proj_mem(traj, top, first_PC, last_PC, aver, evecs):
	flag = False
	print("Projections are calculated while saving memory, it may take some time.")
	ref_aver_str = md.load(aver)
	eigenvecs = load_evecs(evecs)
	mean_vec = ref_aver_str.xyz.astype(np.float64).reshape(1, ref_aver_str.n_atoms * 3)
	if first_PC == None and last_PC == None:
		first_PC = 1
		last_PC = len(eigenvecs)
	elif first_PC == None:
		first_PC = 1
		last_PC = int(last_PC)
	elif last_PC == None:
		first_PC = int(first_PC)
		last_PC = len(eigenvecs)
	else:
		first_PC = int(first_PC)
		last_PC = int(last_PC)
	for frame in md.iterload(traj, top = top, chunk = 100000):
		X = frame.superpose(ref_aver_str).xyz.astype(np.float64).reshape(frame.n_frames, frame.n_atoms * 3) - mean_vec
		if flag == False:
			proj = np.tensordot(X,eigenvecs[first_PC - 1:last_PC],axes = (1,1)).T
			flag = True
		else:
			proj = np.append(proj,np.tensordot(X,eigenvecs[first_PC - 1:last_PC],axes = (1,1)).T,axis = 1)
	return proj, first_PC, last_PC


def main(traj_file, top_file, aver, evecs, memory_flag, first_PC = None, last_PC = None, proj_file = None):
	PATH = os.getcwd() + '/'
	if memory_flag == False:
		proj, first_PC, last_PC = get_proj(load_traj(traj_file, top_file), aver = aver, evecs = evecs, first_PC = first_PC, last_PC = last_PC)
	else:
		proj, first_PC, last_PC = get_proj_mem(traj_file, top_file, aver = aver, evecs = evecs, first_PC = first_PC, last_PC = last_PC)
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
	print('Wrote %s projections in "%s".' % (len(proj),file_out))


# if __name__ == '__main__':
# 	args = sys.argv[1:]
# 	if '-ia' in args and '-f' in args and '-t' in args and '-ievec' in args:
# 		if '-first' in args:
# 			first_PC = int(args[args.index('-first') + 1])
# 		else:
# 			first_PC = None
# 		if '-last' in args:
# 			last_PC = int(args[args.index('-last') + 1])
# 		else:
# 			last_PC = None
# 		if '-op' in args:
# 			proj_file = args[args.index('-op') + 1]
# 		else:
# 			proj_file = None
# 		traj_file = args[args.index('-f') + 1]
# 		top_file = args[args.index('-t') + 1]
# 		aver = args[args.index('-ia') + 1]
# 		evecs = args[args.index('-ievec') + 1]
# 		main(traj_file = traj_file, top_file = top_file, aver = aver, evecs = evecs, first_PC = first_PC, last_PC = last_PC, proj_file = proj_file)
# 	elif '-h' in args or '-help' in args:
# 		print('-f <trajectory file> (file format *.xtc, *trr, etc.)\n-t <topology file> (any file with topology)\n -first <first PC> -last <last PC> \n -ievec <input file with eigenvectors>\n\
# -ia <input file with average structure>\n -op <output file with projections>.')
