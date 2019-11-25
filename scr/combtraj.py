import mdtraj as md
import numpy as np
from multiprocessing import Pool


def load_traj(traj_file, top_file):
	traj = md.load_xtc(traj_file, top_file)
	return traj


def average_structure(traj):
	avg_xyz = traj.xyz.astype(np.float64)
	avg_xyz = avg_xyz.mean(axis = 0, dtype=np.float64)
	avg_traj = md.Trajectory([avg_xyz], traj.top)
	return avg_traj


def main(inp,output_file,conc_file,av_str):
	data = []
	for i in range(0,len(inp),2):
		data.append((inp[i],inp[i+1]))
	with Pool(2) as p:
		trajs = p.starmap(load_traj, data)
	traj = trajs[0]
	length = [len(trajs[i]) for i in range(len(trajs))]
	for i in range(1, len(trajs)):
		traj = traj.join(trajs[i])
		trajs[i] = 0
	traj = traj.superpose(traj[0], parallel = True)
	avg_str = average_structure(traj)
	traj = traj.superpose(avg_str, parallel = True)
	traj[:length[0]].save_xtc(conc_file[:conc_file.rfind('.')]+str(1)+'.xtc')
	average_structure(traj[:length[0]]).save(av_str[:av_str.rfind('.')]+str(1)+'.pdb')
	start_ = length[0]
	end_ = length[0]
	for i in range(1, len(length)):
		traj[start_:end_ + length[i]].save_xtc(\
			conc_file[:conc_file.rfind('.')]+str(i+1)+'.xtc')
		average_structure(traj[start_:end_ + length[i]]).save(\
			av_str[:av_str.rfind('.')]+str(i+1)+'.pdb')
		start_ += length[i]
		end_ += length[i]
	avg_str = average_structure(traj)
	traj.save_xtc(output_file)
	avg_str.save(av_str)

