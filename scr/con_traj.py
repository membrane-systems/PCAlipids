import mdtraj as md
import numpy as np
from multiprocessing import Pool
import sys
import os
import time
import gc
#from numba import jit


def load_traj(traj_file, traj_top, lipid_resname, stride, sf, ef, sel = "", max_frames = None):
	# create topology
	topol = md.load(traj_top).topology
	# select atoms for loading
	if sel == "":
		ailist = topol.select('not water and not type W H Hs WT4 NaW KW CLW MgW and resname %s' \
			% lipid_resname)
	else:
		ailist = topol.select(sel)
	# load trajectory
	return md.load(traj_file, top = traj_top, stride = stride, \
		atom_indices = ailist)[sf:ef] if ef != -1 else \
                md.load(traj_file, top = traj_top, stride = stride, \
                atom_indices = ailist)[sf:]

def concat(traj, lipid):
	return traj.atom_slice(traj.topology.select('resid %s' % lipid))


def average_structure(traj):
	avg_xyz = traj.xyz.mean(axis = 0, dtype=np.float64)
	avg_traj = md.Trajectory([avg_xyz], traj.top)
	return avg_traj

def main(file_1, file_2, stride, sf, ef, \
	out_traj, out_top, file_3, lipid_resname, lipid_sel):
	
	if not file_1 or not file_2 or not lipid_resname:
		print("Trajectory and topology files have to be provided.\n\
Name of residue of interest have to be indicated.\n\
Run pcalipids.py concat -h for help")
		return 0

	# load trajectory
	PATH = os.getcwd() + '/'
	traj = load_traj(PATH + file_1, PATH + file_2, \
		lipid_resname = lipid_resname, stride = stride, sf = sf, ef = ef, sel = lipid_sel)
	print("Loaded trajectory "+file_1)

	# get number of frames, lipids and atoms
	nframes = traj.n_frames
	nlip    = len(list(traj.topology.residues))
	na      = traj.xyz.shape[1]//nlip
	
	# create concatenated trajectory
	trajxyz = traj.xyz
	trajxyz = trajxyz.reshape((nframes,nlip,na,3))
	trajxyz = np.swapaxes(trajxyz,0,1)
	trajxyz = trajxyz.reshape((nframes*nlip,na,3))

	traj = md.Trajectory(trajxyz,traj[0].atom_slice(range(na)).topology)
	
	print("Created concatenated trajectory")

	# aling frames in the concatenated trajectory
	if file_3:
		# if reference structure is provided:
		# align all frames to the reference structure

		ref_traj = md.load(file_3)
		traj.superpose(ref_traj, parallel = True)
	else:
		# if reference frame is not provided
		# align all frames to the first frame of the first lipid
		# then align all frames to the average structure of the aligned traj
		ref_traj = traj[0]
		traj.superpose(ref_traj, parallel = True)
		avg_str = average_structure(traj)
		traj.superpose(avg_str, parallel = True)
	print('The trajectory was aligned, reference = %s.' % file_3) if file_3 else \
	print('The trajectory was aligned to the average structure')
	
	# save concatenated trajectory
	traj.save(PATH + out_traj)
	print('Concatenated trajectory was created, saved in %s' % out_traj)
	
	# save average structure
	avg_str = average_structure(traj)
	avg_str.save(out_top)
	print('Average structure saved in "%s"' % out_top)
	return
