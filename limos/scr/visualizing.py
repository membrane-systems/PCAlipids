import matplotlib.pyplot as pyplot
import mdtraj as md
import numpy as np
import sys


def extreme_proj(proj_file, PC):
	with open(proj_file, 'r') as file:
		line = file.readline()
		i = 0
		while i < int(PC) - 1:
			line = file.readline()
			if line.find('&') != -1:
				i += 1

		min_ = 10 ** 10
		max_ = (-1) * 10 **10
		line = file.readline()
		while line.find('&') == -1:
			if line.find('@') == -1:
				print(line)
				proj_val = float(line.split()[1])
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
			line = file.realdine()
			i += 1
		eigvec = []
		for i in range(len(line.split())):
			eigvec.append(np.float64(line.split()[i]))

	return np.array(eigvec, dtype = np.float64).T


def create_extreme(eigvec, avg_str_xyz, max_, min_, avg_str):
	extr_max_xyz = [eigvec[i] * max_ + avg_str_xyz[i] for i in range(len(eigvec))]
	extr_max_xyz = np.array(extr_max_xyz).reshape(len(extr_max_xyz) // 3, 3)
	print(extr_max_xyz)
	extr_min_xyz = [eigvec[i] * min_ + avg_str_xyz[i] for i in range(len(eigvec))]
	extr_min_xyz = np.array(extr_min_xyz).reshape(len(extr_min_xyz) // 3, 3)
	print(extr_min_xyz)
	extr_max_traj = md.Trajectory([extr_max_xyz], avg_str.top)
	extr_min_traj = md.Trajectory([extr_min_xyz], avg_str.top)
	extr_max_traj.save('extreme_1_max.pdb')
	extr_min_traj.save('extreme_1_min.pdb')


def interpolate(extreme_file_max, extreme_file_min):
	extr_max_traj = md.load(extreme_file_max)
	extr_min_traj = md.load(extreme_file_min)

	extr_max_xyz = extr_max_traj.xyz.astype(np.float64).reshape(1,extr_max_traj.n_atoms*3).T
	extr_min_xyz = extr_min_traj.xyz.astype(np.float64).reshape(1,extr_min_traj.n_atoms*3).T
	print(extr_max_xyz - extr_min_xyz)
	delta_vec = np.array((extr_max_xyz - extr_min_xyz), dtype = np.float64)
	print(delta_vec)
	inter_str = extr_min_traj
	for i in range(1, 20):
		inter_str = inter_str.join(md.Trajectory([(extr_min_xyz + delta_vec * i * 0.05).reshape(len(extr_min_xyz) // 3, 3)], extr_max_traj.top))
	inter_str = inter_str.join(extr_max_traj)
	inter_str.save('extreme1.pdb')


def main(proj_file, PC, pdb_file, eigvec_file):
	max_, min_ = extreme_proj(proj_file, PC)
	print(max_, min_)
	avg_str_xyz, avg_str = get_average_struct(pdb_file)
	eigvec = get_eigvec(eigvec_file, PC)
	print(eigvec)
	create_extreme(eigvec, avg_str_xyz, max_, min_, avg_str)
	interpolate('extreme_1_max.pdb', 'extreme_1_min.pdb')


if __name__ == '__main__':
	args = sys.argv[1:]
	main(args[args.index('-p') + 1], args[args.index('-n') + 1], args[args.index('-pdb') + 1], args[args.index('-e') + 1])