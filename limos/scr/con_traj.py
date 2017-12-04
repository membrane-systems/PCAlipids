import mdtraj as md
import numpy as np
from multiprocessing import Pool
import sys
import os


def load_traj(traj_file, traj_top, lipid_resname, max_frames = None):
	print('Loading trajctory..')
	traj = md.load(traj_file, top = traj_top)
	print('Removing solvent..')
	if lipid_resname == None:
		traj = traj.atom_slice(traj.topology.select('not water and not type W H Hs'))
	else:
		traj = traj.atom_slice(traj.topology.select('not water and not type W H Hs and resname %s' % lipid_resname))
	return traj


def concat(traj, lipid):
	sin_lip_traj = traj.atom_slice(traj.topology.select('resid %s' % lipid))
	return sin_lip_traj


def average_structure(traj):
	avg_xyz = traj.xyz.astype(np.float64)
	avg_xyz = avg_xyz.mean(axis = 0, dtype=np.float64)
	avg_traj = md.Trajectory([avg_xyz], traj.top)
	return avg_traj


def main(file_1, file_2, file_3 = None, lipid_resname = None):
	PATH = os.getcwd() + '/'
	traj = load_traj(PATH + file_1, PATH + file_2, lipid_resname = lipid_resname)
	N = traj[0].n_residues
	#lipids = [(traj,i)  for i in range(N)]
	#with Pool(2) as p:
	#	trajs_ = p.starmap(concat, lipids)
	trajs_ = []
	for i in range(N):
		trajs_.append(concat(traj, i))	
	traj = trajs_[0]
	for i in range(1, N):
		traj = traj.join(trajs_[i])
		trajs_[i] = 0
	if file_3 != None:
		ref_traj = md.load(file_3)
	else:
		ref_traj = traj[0]
		print('The trajectory was aligned, reference = first frame.')
	traj = traj.superpose(ref_traj)
	traj.save_xtc(PATH + 'concatenated.xtc')
	print('Concatenated trajectory was created, saved in "concatenated.xtc"')
	avg_str = average_structure(traj)
	avg_str.save('average.pdb')
	print('Average structure saved in "average.pdb"')
	return 'concatenated.xtc', 'average.pdb'


if __name__ == '__main__':
	args = sys.argv[1:]
	if '-f' in args and '-t' in args and '-r' in args and '-l' in args:
		main(args[args.index('-f') + 1], args[args.index('-t') + 1], args[args.index('-r') + 1], args[args.index('-l') + 1])
	elif '-f' in args and '-t' in args and '-r' not in args:
		print('No reference file supplied. The first frame of trajectory will be used for alignment.')
		main(args[args.index('-f') + 1], args[args.index('-t') + 1])
	elif '-f' in args and '-t' in args and '-r' in args and '-l' not in args:
		main(args[args.index('-f') + 1], args[args.index('-t') + 1], args[args.index('-r') + 1])
	elif '-f' in args and '-t' in args and '-r' not in args and '-l' not in args:
		main(args[args.index('-f') + 1], args[args.index('-t') + 1])
	elif '-h' not in args:
		print('Missing parameters, try -h for flags\n')
	else:
		print('-f <trajectory file> (file format *.xtc, *trr)\n-t <topology file> (any file with topology)\n-r <reference traj file> (any topology file). If not supplied, \
			the first frame of trajectory will be used for alignment\n -l <lipid type> (example: -l DPPC).')

