import matplotlib.pyplot as pyplot
import mdtraj as md
import numpy as np
import sys


def extreme_proj(proj_file):
	with open(proj_file, 'r') as file:
		min_ = 10 ** 10
		max_ = (-1) * 10 **10

		line = file.readline()

		while line.find('&') == -1:
			if line.find('@') > -1:
				line = file.readline()
				continue
			proj_val = float(line)
			max_ = max(max_, proj_val)
			min_ = min(min_, proj_val)
			line = file.readline()
	return max_, min_


def get_average_struct(pdb_file):
	avg_str = md.load(pdb_file)
	avg_str_xyz = avg_str.xyz.astype(np.float64).reshape(1,avg_str.n_atoms*3).T
	buff = [avg_str_xyz[i] for i in range(len(avg_str_xyz))]
	buff = np.array(buff, dtype = np.float64)
	return buff, avg_str


def get_eigvec(eigvec_file, PC):
	with open(eigvec_file, 'r') as file:
		i = 1
		line = file.readline()
		while i < int(PC):
			line = file.readline()
			i += 1
		eigvec = []
		for i in range(len(line.split())):
			eigvec.append(np.float64(line.split()[i]))

	return np.array(eigvec, dtype = np.float64).T


def create_extreme(eigvec, avg_str_xyz, max_, min_, avg_str, PC):
	extr_max_xyz = [eigvec[i] * max_ + avg_str_xyz[i] for i in range(len(eigvec))]
	extr_max_xyz = np.array(extr_max_xyz).reshape(len(extr_max_xyz),1)
	
	extr_min_xyz = [eigvec[i] * min_ + avg_str_xyz[i] for i in range(len(eigvec))]
	extr_min_xyz = np.array(extr_min_xyz).reshape(len(extr_min_xyz),1)
	
	delta_vec = np.array((extr_max_xyz - extr_min_xyz), dtype = np.float64)

	extr_max_traj = md.Trajectory([extr_max_xyz.reshape(len(extr_max_xyz)//3,3)], \
		avg_str.top)
	extr_min_traj = md.Trajectory([extr_min_xyz.reshape(len(extr_min_xyz)//3,3)], \
		avg_str.top)

	inter_str = extr_min_traj
	for i in range(1, 21):
		inter_str = inter_str.join(md.Trajectory([\
			(extr_min_xyz + delta_vec * i * 0.05).reshape(len(extr_min_xyz) // 3, 3)], \
			extr_max_traj.top))

	inter_str.save('extreme_'+str(PC)+'.pdb')
	print("Projection of "+str(PC)+" principal component is saved to extreme_"+str(PC)+".pdb")


def main(proj_file, pdb_file, eigvec_file):
	# Get PC #
	PC = int(proj_file[proj_file.rfind('_')+1:proj_file.rfind('.')])

	# Find maximal and minimal projection values
	max_, min_ = extreme_proj(proj_file)
	
	# Read average structure
	avg_str_xyz, avg_str = get_average_struct(pdb_file)

	# Get eigenvector for PC
	eigvec = get_eigvec(eigvec_file, PC)
	
	# Create extreme projections
	create_extreme(eigvec, avg_str_xyz, max_, min_, avg_str, PC)
#	interpolate('extreme_1_max.pdb', 'extreme_1_min.pdb')

