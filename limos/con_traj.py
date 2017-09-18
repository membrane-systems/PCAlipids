import mdtraj as md
from multiprocessing import Pool
import sys
import os


def load_traj(traj_file, traj_top, max_frames = None):
	print('Loading trajctory..')
	traj = md.load_xtc(traj_file, top = traj_top)
	print('Removing solvent..')
	traj = traj.atom_slice(traj.topology.select('not water and not type W H Hs'))
	return traj


def concat(traj, lipid):
	sin_lip_traj = traj.atom_slice(traj.topology.select('resid %s' % lipid))
	return sin_lip_traj


def main(file_1, file_2, file_3 = None):
	PATH = os.getcwd() + '/'
	traj = load_traj(PATH + file_1, PATH + file_2)
	N = traj[0].n_residues
	lipids = [(traj,i)  for i in range(N)]
	with Pool(2) as p:
		trajs_ = p.starmap(concat, lipids)
	traj = trajs_[0]
	for i in range(1, N):
		traj = traj.join(trajs_[i])
		trajs_[i] = 0
	print('Concatenated trajectory was created.')
	if file_3 != None:
		ref_traj = md.load(PATH + file_3)
	else:
		ref_traj = traj[0]
	traj = traj.superpose(ref_traj)
	print('The trajectory was aligned.')
	traj.save_xtc(PATH + 'full.xtc')


if __name__ == '__main__':
	args = sys.argv[1:]
	if '-f' in args and '-t' in args and '-r' in args:
		main(args[args.index('-f') + 1], args[args.index('-t') + 1], args[args.index('-r') + 1])
	elif '-f' in args and '-t' in args and '-r' not in args:
		print('No reference file supplied. The first frame of trajectory will be used for alignment.')
		main(args[args.index('-f') + 1], args[args.index('-t') + 1])
	elif '-h' not in args:
		print('Missing parameters, try -h for flags\n')
	else:
		print('-f <trajectory file> (file format *.xtc)\n-t <topology file> (any file with topology)\n-r <reference traj file> (any file with 1 ref frame). If not supplied, \
			the first frame of trajectory will be used for alignment.')
