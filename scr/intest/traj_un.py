import mdtraj as md
import numpy as np
import sys
import os
from multiprocessing import Pool
import time


def load_traj(traj_file, top_file):
	traj = md.load_xtc(traj_file, top_file)
	return traj


def main(data):
	with Pool(2) as p:
		trajs = p.starmap(load_traj, data)
	traj = trajs[0]
	for i in range(1, len(trajs)):
		traj.join(trajs[i])
		trajs[i] = 0
	traj.save_xtc('associated.xtc')


if __name__ == '__main__':
	args = sys.argv[1:]
	i = 0
	data = []
	while i < len(args):
		if args[i] == '-fs':
			data.append((args[i + 1], args[i + 2]))
			i += 3
	main(data)