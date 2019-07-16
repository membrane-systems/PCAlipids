import mdtraj as md
import numpy as np
from scipy.stats import skew
import sys
import os


def load_traj(traj_file, top_file):
	traj = md.load(traj_file, top = top_file)
	return traj


def load_evecs(evecs):
	eigenvecs = []
	file = open(evecs, 'r')
	for line in file:
		a = line.split()
		for i in range(len(a)):
			a[i] = np.float64(a[i])
		a = np.array(a, dtype = np.float64)
		eigenvecs.append(a)
	file.close()
	eigenvecs = np.array(eigenvecs, dtype = np.float64)
	return eigenvecs


def get_proj(traj, first_PC, last_PC, aver, evecs, filename):
	ref_aver_str = md.load(aver)
	traj = traj.superpose(ref_aver_str)
	mean_vec = ref_aver_str.xyz.astype(np.float64)
	x_std = traj.xyz.astype(np.float64).reshape(traj.n_frames,traj.n_atoms*3).T
	x_std = x_std - np.array([mean_vec.reshape(traj.n_atoms*3, )]).T
	eigenvecs = load_evecs(evecs)
	if first_PC == None and last_PC == None:
		first_PC = 1
		last_PC = len(eigenvecs)
	elif first_PC == None:
		first_PC = 1
		last_PC = int(last_PC)
	elif last_PC == None:
		first_PC = int(first_PC)
		last_PC = len(eigenvecs)
	else:
		first_PC = int(first_PC)
		last_PC = int(last_PC)
	files = []
	for i in range(first_PC,last_PC + 1):
		file = open(filename[:filename.rfind('.')]+'_'+str(i)+filename[filename.rfind('.'):],'w')
		file.write('@    title "Projection %s"\n' % (i))
		files.append(file)
	for i in range(first_PC - 1, last_PC):
		proj = (x_std).T.dot(eigenvecs[i])
		files[i].write(''.join(('     ' + str(proj[j]) + '\n') for j in range(len(proj))))
	for i in range(len(files)):
		files[i].write('&\n')
		files[i].close
	print('Wrote %s projections in "%s".' % (len(files), files[0].name[:files[0].name.rfind('_')] + "*"))

def get_proj_mem(traj, top, first_PC, last_PC, aver, evecs, filename):
	print("Projections are calculated while saving memory, it may take some time.")
	
	ref_aver_str = md.load(aver)
	na = ref_aver_str.n_atoms
	eigenvecs = load_evecs(evecs)
	mean_vec = ref_aver_str.xyz.astype(np.float64).reshape(1, ref_aver_str.n_atoms * 3)
	
	first_PC = first_PC
	last_PC = last_PC
	
	files = []
	
	chPC1 = (first_PC == 1)
	
	# if chPC1 - we need to orient it correctly

	if chPC1:
		proj0 = [] #create the array for pc1 projections

	# prepare files
	for i in range(first_PC,last_PC + 1):
		file = open(filename[:filename.rfind('.')]+'_'+str(i)+filename[filename.rfind('.'):],'w')
		file.write('@    title "Projection %s"\n' % (i))
		files.append(file)
	
	for frame in md.iterload(traj, top = top, chunk = 100000): 
		# superpose trajectory
		X = frame.superpose(ref_aver_str).xyz.astype(np.float64)\
		.reshape(frame.n_frames, frame.n_atoms * 3) - mean_vec
	
		# calculate projections for chunk
		proj = np.tensordot(X,eigenvecs[first_PC - 1:last_PC],axes = (1,1)).T 
		
		for i in range(len(files)):
			# Find wether the distribution is skewed to the correct side
			if chPC1 and i == 0:
				proj0.append(proj[i]) #save projection for PC1 into special array
			else:
				files[i].write(''.join(('     ' + str(proj[i][j]/(na**0.5)) + '\n') \
				for j in range(len(proj[i])))) # write projection for other PCs directly to file
	# Write projection on PC1 to file		
	if chPC1:
		proj0=np.concatenate(np.array(proj0),axis=None)
		k = skew(proj0)/abs(skew(proj0)) # calculate the correct orientation

		# if k < 0 we then need to rewright eigenvector
		# it's direction has to be changed
		if k < 0:
			eigenvecs[0] = -eigenvecs[0]
			np.savetxt(evecs,eigenvecs,fmt='%20.17g')
			print('Wrote eigenvectors in "%s"' % evecs)

		files[0].write(''.join(('     ' + str(k*proj0[j]/(na**0.5)) + '\n') \
				for j in range(len(proj0))))
	
	# Correct endings to files	
	for i in range(len(files)):
		files[i].write('&\n')
		files[i].close
	
	print('Wrote %s projections in "%s".' % (len(files), files[0].name[:files[0].name.rfind('_')] + "*"))


def main(traj_file, top_file, aver, evecs, first_PC, last_PC, proj_file):
	PATH = os.getcwd() + '/'
	get_proj_mem(traj_file, top_file, aver = aver, evecs = evecs, \
		first_PC = first_PC, last_PC = last_PC, filename = PATH + proj_file)
