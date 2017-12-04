import mdtraj as md
import numpy as np
from multiprocessing import Pool
import sys
import os


def load_traj(traj_file, traj_top, lipid_resname, stride, sf, max_frames = None):
	print('Loading trajectory..')
	if sf == None:
		traj = md.load(traj_file, top = traj_top, stride = stride)
	else:
		traj = md.load(traj_file, top = traj_top, stride = stride)[sf:]
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


def main(file_1, file_2, stride, sf, out_traj, out_top, file_3 = None, lipid_resname = None):
	PATH = os.getcwd() + '/'
	traj = load_traj(PATH + file_1, PATH + file_2, lipid_resname = lipid_resname, stride = stride, sf = sf)
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
	if out_traj == None:
		traj.save(PATH + 'concatenated.xtc')
		print('Concatenated trajectory was created, saved in "concatenated.xtc"')
	else:
		traj.save(PATH + out_traj)
		print('Concatenated trajectory was created, saved in "%s"' % out_traj)
	avg_str = average_structure(traj)
	if out_top == None:
		avg_str.save('average.pdb')
		print('Average structure saved in "average.pdb"')
	else:
		avg_str.save(out_top)
		print('Average structure saved in "%s"' % out_top)
	return


if __name__ == '__main__':
	args = sys.argv[1:]
	if '-stride' in args:
		stride = int(args[args.index('-stride') + 1])
	else:
		stride = None
	if '-sf' in args:
		sf = int(args[args.index('-sf') + 1])
	else:
		sf = None
	if '-oc' in args:
		out_traj = args[args.index('-oc') + 1]
	else:
		out_traj = None
	if '-oa' in args:
		out_top = args[args.index('-oa') + 1]
	else:
		out_top = None
	if '-f' in args and '-t' in args and '-r' in args and '-l' in args:
		main(args[args.index('-f') + 1], args[args.index('-t') + 1], args[args.index('-r') + 1], args[args.index('-l') + 1], stride = stride, sf = sf, out_traj = out_traj, out_top = out_top)
	elif '-f' in args and '-t' in args and '-r' not in args:
		print('No reference file supplied. The first frame of trajectory will be used for alignment.')
		main(args[args.index('-f') + 1], args[args.index('-t') + 1], stride = stride, sf = sf, out_traj = out_traj, out_top = out_top)
	elif '-f' in args and '-t' in args and '-r' in args and '-l' not in args:
		main(args[args.index('-f') + 1], args[args.index('-t') + 1], args[args.index('-r') + 1], stride = stride, sf = sf, out_traj = out_traj, out_top = out_top)
	elif '-f' in args and '-t' in args and '-r' not in args and '-l' not in args:
		main(args[args.index('-f') + 1], args[args.index('-t') + 1], stride = stride, sf = sf, out_traj = out_traj, out_top = out_top)
	elif '-h' not in args:
		print('Missing parameters, try -h for flags\n')
	else:
		print('-f <trajectory file> (file format *.xtc, *trr)\n-t <topology file> (any file with topology)\n-r <reference traj file> (any topology file). If not supplied, \
			the first frame of trajectory will be used for alignment\n -l <lipid type> (example: -l DPPC)\n -stride <positive integer; step of reading frames>\n \
			 -sf <time in ps; number to determine from which frame to read the trajectory>\n -oc <output trajectory file>\n -oa <output topology file>\n')

