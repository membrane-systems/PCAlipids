import mdtraj as md
import numpy as np
import sys


def lipic(traj_file, top_file, stride, toAl):
	# load trajectory
	traj = md.load_xtc(traj_file, top = top_file, stride = stride)
	# align if needed
	if toAl:
		# align to the first frame
		traj.superpose(traj[0])

		# get average structure
		avg_xyz = traj.xyz.astype(np.float64).\
		mean(axis = 0, dtype=np.float64)
		avg_traj = md.Trajectory([avg_xyz], traj.top)

		# second step alignment to the average structure
		traj.superpose(avg_traj[0])
		
	# calculate average structure
	avg_xyz = traj.xyz.astype(np.float64).\
	mean(axis = 0, dtype=np.float64)
	avg_traj = md.Trajectory([avg_xyz], traj.top)

	# add average structure to the trajectory
	traj = traj.join(avg_traj)
	
	return traj


def main(traj_file, top_file, stride, conf_out, toAl):
	if not traj_file or not top_file:
		print("Trajectory and topology have to be provided. \n\
			Please run pcalipids.py conspace -h for help")
		return
	# create pdb file with conformations
	traj = lipic(traj_file, top_file, stride, toAl)
	traj.save_pdb(conf_out)
	print('Conformations saved in file "%s"' % conf_out)
