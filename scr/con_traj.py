import mdtraj as md
import numpy as np
from multiprocessing import Pool
import sys
import os
#from numba import jit


def load_traj(traj_file, traj_top, lipid_resname, stride, sf, ef = None, max_frames = None):
	print('Loading trajectory..')
	if sf == None and ef == None:
		traj = md.load(traj_file, top = traj_top, stride = stride)
	else:
		traj = md.load(traj_file, top = traj_top, stride = stride)[sf:ef]
	print('Removing solvent..')
	if lipid_resname == None:
		traj = traj.remove_solvent()
		traj = traj.atom_slice(traj.topology.select('not water and not type W H Hs WT4 NaW KW CLW MgW'))
	else:
		traj = traj.remove_solvent()
		traj = traj.atom_slice(traj.topology.select('not water and not type W H Hs WT4 NaW KW CLW MgW and resname %s' % lipid_resname))
	print(traj)
	# # table, bonds  = traj.topology.to_dataframe()
	# table.loc[:, ('chainID')] = 0
	# table.loc[:, ('segmentID')] = 'A'
	# topology = md.Topology.from_dataframe(table, bonds)
	# new_traj = md.Trajectory(traj.xyz.astype(np.float64), topology = topology)
	# traj = 0
	return traj


def concat(traj, lipid):
	return traj.atom_slice(traj.topology.select('resid %s' % lipid))


def average_structure(traj):
	avg_xyz = traj.xyz.astype(np.float64)
	avg_xyz = avg_xyz.mean(axis = 0, dtype=np.float64)
	avg_traj = md.Trajectory([avg_xyz], traj.top)
	return avg_traj


#@jit
def main(file_1, file_2, stride, sf, ef, out_traj, out_top, file_3 = None, lipid_resname = None):

	PATH = os.getcwd() + '/'
	traj = load_traj(PATH + file_1, PATH + file_2, lipid_resname = lipid_resname, stride = stride, sf = sf, ef = ef)
	N = traj[0].n_residues
	trajs_ = []
	i = 0
	j = 0
	while j < N:
		try:
			trajs_.append(concat(traj, i))	
			i+= 1
			j+= 1
		except IndexError:
			i += 1
	traj = trajs_[0]
	for i in range(1, N):
		traj = traj.join(trajs_[i])
		trajs_[i] = 0
	if file_3 != None:
		ref_traj = md.load(file_3)
		print('The trajectory was aligned, reference = %s.' % file_3)
	else:
		ref_traj = traj[0]
		print('The trajectory was aligned, reference = first frame.')
	traj = traj.superpose(ref_traj, parallel = True)
	if out_traj == None:
		traj.save(PATH + 'concatenated.xtc')
		print('Concatenated trajectory was created, saved in "concatenated.xtc"')
	else:
		traj.save(PATH + out_traj)
		print('Concatenated trajectory was created, saved in "%s"' % out_traj)
	avg_str = average_structure(traj)
	traj = traj.superpose(avg_str, parallel = True)
	avg_str = average_structure(traj)
	if out_top == None:
		avg_str.save('average.pdb')
		print('Average structure saved in "average.pdb"')
	else:
		avg_str.save(out_top)
		print('Average structure saved in "%s"' % out_top)
	return
