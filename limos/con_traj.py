import time
import mdtraj as md
from multiprocessing import Pool
import sys
import os


def load_traj(traj_file, traj_top, max_frames = None):
	traj = md.load_xtc(traj_file, top = traj_top)
	traj = traj.atom_slice(traj.topology.select('not water and not type W H Hs'))
	return traj


def concat(traj, lipid):
	sin_lip_traj = traj.atom_slice(traj.topology.select('resid %s' % lipid))
	return sin_lip_traj


def main(file_1, file_2):
	PATH = os.getcwd()
	traj = load_traj(PATH + '/' + file_1, PATH + '/' + file_2)
	N = traj[0].n_residues
	lipids = [(traj,i)  for i in range(N)]
	with Pool(2) as p:
		trajs_ = p.starmap(concat, lipids)
	traj = trajs_[0]
	for i in range(1, N):
		traj = traj.join(trajs_[i])
		trajs_[i] = 0
	traj.save_xtc(PATH + '/full.xtc')
	

if __name__ == '__main__':
	args = sys.argv[1:]
	if '-f' in args and '-t' in args:
		main(args[args.index('-f') + 1], args[args.index('-t') + 1])
	elif '-h' not in args:
		print('Missing parameters, try -h for flags\n')
	else:
		print('-f <trajectory file> (file format *.xtc)\n-t <topology file> (any file with topology)')