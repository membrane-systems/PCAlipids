#!/usr/bin/env python3

import PCA_old
import con_traj
import sys
import kss_time_scale
import autocorr
import eigenvecdot
import lipic
import proj_dist_s
import pearson
import proj_parse
import visualizing
import PCA
import project
import projdistm
import eigenvalues
import combtraj
import reltime
import timescalespic

def main(args):
	if  args[0] == 'concat':
		main = con_traj.main
		if '-stride' in args:
			stride = int(args[args.index('-stride') + 1])
		else:
			stride = None
		if '-sf' in args:
			sf = int(args[args.index('-sf') + 1])
		else:
			sf = None
		if '-ef' in args:
			ef = int(args[args.index('-ef') + 1])
		else:
			ef = None
		if '-oc' in args:
			out_traj = args[args.index('-oc') + 1]
		else:
			out_traj = None
		if '-oa' in args:
			out_top = args[args.index('-oa') + 1]
		else:
			out_top = None
		if '-f' in args and '-t' in args and '-r' in args and '-l' in args:
			main(args[args.index('-f') + 1], args[args.index('-t') + 1], file_3 = args[args.index('-r') + 1], lipid_resname = args[args.index('-l') + 1], stride = stride, sf = sf, ef = ef, out_traj = out_traj, out_top = out_top)
		
		elif '-f' in args and '-t' in args and '-l' in args and '-r' not in args:
			print('No reference file supplied. The first frame of trajectory will be used for alignment.')
			main(args[args.index('-f') + 1], args[args.index('-t') + 1], lipid_resname = args[args.index('-l') + 1], stride = stride, sf = sf, ef = ef, out_traj = out_traj, out_top = out_top)
		
		# elif '-f' in args and '-t' in args and '-r' in args and '-l' not in args:
		# 	main(args[args.index('-f') + 1], args[args.index('-t') + 1], args[args.index('-r') + 1], stride = stride, sf = sf, out_traj = out_traj, out_top = out_top)
		
		# elif '-f' in args and '-t' in args and '-r' not in args and '-l' not in args:
		# 	main(args[args.index('-f') + 1], args[args.index('-t') + 1], stride = stride, sf = sf, out_traj = out_traj, out_top = out_top)
		
		elif '-h' not in args:
			print('Missing parameters, try -h for flags\n')
		
		else:
			print(' -f <trajectory file> (file format *.xtc, *trr)\n -t <topology file> (any file with topology)\n -r <reference traj file> (any topology file). If not supplied, \
the first frame of trajectory will be used for alignment\n -l <lipid type> (example: -l DPPC)\n -stride <positive integer; step of reading frames>\n \
-sf <time in ps; number to determine from which frame to read the trajectory>\n -ef <time in ps; number to determine to which frame to read the trajectory>\n -oc <output trajectory file>\n -oa <output topology file>\n')

	
	elif args[0] == 'covar':
		main = PCA.main
		if '-oeval' in args:
			val_file = args[args.index('-oeval') + 1]
		else:
			val_file = None
		if '-oevec' in args:
			vec_file = args[args.index('-oevec') + 1]
		else:
			vec_file = None
		if '-ocov' in args:
			cov_file = args[args.index('-ocov') + 1]
		else:
			cov_file = None
		if '-invertPC1' in args:
			invert = args[args.index('-invertPC1') + 1]
			if invert == 'False' or invert == '0':
				invert = False
			else:
				invert = True
		else:
			invert = False
		if '-emem' in args:
			memory_flag = True
		else:
			memory_flag = False
		if '-f' in args and '-t' in args:
			main(args[args.index('-f') + 1], args[args.index('-t') + 1], val_file = val_file, vec_file = vec_file, cov_file = cov_file, invert = invert, memory_flag = memory_flag)
		elif '-h' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-f <trajectory file> (file format *.xtc)\n-t <topology file> (any file with topology)\n -first <first PC> -last <last PC> \n -oeval <output file with eigenvalues>\n\
-oevec <output file with eigenvectors>\n -ocov <output file with covariance matrix>\n -invertPC1 <invert PC1>\n -emem <memory economy>')


	elif args[0] == 'project':
		main = project.main
		if '-ia' in args and '-f' in args and '-t' in args and '-ievec' in args:
			if '-first' in args:
				first_PC = int(args[args.index('-first') + 1])
			else:
				first_PC = None
			if '-last' in args:
				last_PC = int(args[args.index('-last') + 1])
			else:
				last_PC = None
			if '-op' in args:
				proj_file = args[args.index('-op') + 1]
			else:
				proj_file = None
			if '-emem' in args:
				memory_flag = True
			else:
				memory_flag = False
			traj_file = args[args.index('-f') + 1]
			top_file = args[args.index('-t') + 1]
			aver = args[args.index('-ia') + 1]
			evecs = args[args.index('-ievec') + 1]
			main(traj_file = traj_file, top_file = top_file, aver = aver, evecs = evecs, memory_flag = memory_flag, first_PC = first_PC, last_PC = last_PC, proj_file = proj_file)
		elif '-h' in args or '-help' in args:
			print('-f <trajectory file> (file format *.xtc, *trr, etc.)\n-t <topology file> (any file with topology)\n -first <first PC> -last <last PC> \n -ievec <input file with eigenvectors>\n\
-ia <input file with average structure>\n -op <output file with projections>\n -emem <memory economy>.')
		else:
			print('Missing parameters, try -h for flags\n')

	
	elif args[0] == 'projdist':
		main = proj_dist_s.main
		if '-p' in args:
			filename = args[args.index('-p') + 1]
			PC = int(filename.split('_')[1].split('.')[0])
			main(filename, PC)


		elif '-pr' in args and '-pr' in args:
			files = args[args.index('-pr') + 1]
			file_start = files[:files.find('-')]
			file_end = files[files.find('-') + 1:]
			start = int(file_start[file_start.rfind('_') + 1:file_start.rfind('.')])
			end = int(file_end[file_end.rfind('_') + 1:file_end.rfind('.')])
			file_mask = file_start[:file_start.rfind('_')] + file_start[file_start.rfind('.'):]
			filenames = [file_mask[:file_mask.find('.')] + "_" + str(i) + file_mask[file_mask.find('.'):] for i in range(start, end + 1)]

			PCs = [int(filenames[i].split('_')[1].split('.')[0]) for i in range(len(filenames))]
			for i in range(len(PCs)):
				main(filenames[i],PCs[i])

		elif '-h' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-p <projection file> (file format *.xvg)\n-npc - number of component')

	
	elif args[0] == 'ksst':
		main = kss_time_scale.main
		if '-p' in args and '-ln' in args and '-dt' in args and '-o' in args and '-pr' not in args:
			start = args.index('-p')
			end = min(args.index('-o'), args.index('-ln'))
			filenames = [args[i] for i in range(start + 1, end)]
			main(filenames, int(args[args.index('-ln') + 1]), float(args[args.index('-dt') + 1]), args[args.index('-o') + 1])
		elif '-p' in args and '-ln' in args and '-dt' in args and '-o' not in args and '-pr' not in args:
			# print('No output file supplied. Data will be written in "KSS_relaxtime_vs_PC.xvg"')
			filenames = [args[i] for i in range(args.index('-p') + 1, min(args.index('-dt'), args.index('-ln')))]
			main(filenames, int(args[args.index('-ln') + 1]), float(args[args.index('-dt') + 1]))
		elif '-pr' in args and '-ln' in args and '-dt' in args and '-o' in args and '-p' not in args:
			files = args[args.index('-pr') + 1]
			file_start = files[:files.find('-')]
			file_end = files[files.find('-') + 1:]
			start = int(file_start[file_start.rfind('_') + 1:file_start.rfind('.')])
			end = int(file_end[file_end.rfind('_') + 1:file_end.rfind('.')])
			file_mask = file_start[:file_start.rfind('_')] + file_start[file_start.rfind('.'):]
			filenames = [file_mask[:file_mask.find('.')] + "_" + str(i) + file_mask[file_mask.find('.'):] for i in range(start, end + 1)]
			main(filenames, int(args[args.index('-ln') + 1]), float(args[args.index('-dt') + 1]), args[args.index('-o') + 1])
		elif '-pr' in args and '-ln' in args and '-dt' in args and '-o' not in args and '-p' not in args: 
			files = args[args.index('-pr') + 1]
			file_start = files[:files.find('-')]
			file_end = files[files.find('-') + 1:]
			start = int(file_start[file_start.rfind('_') + 1:file_start.rfind('.')])
			end = int(file_end[file_end.rfind('_') + 1:file_end.rfind('.')])
			file_mask = file_start[:file_start.rfind('_')] + file_start[file_start.rfind('.'):]
			filenames = [file_mask[:file_mask.find('.')] + "_" + str(i) + file_mask[file_mask.find('.'):] for i in range(start, end + 1)]
			main(filenames, int(args[args.index('-ln') + 1]), float(args[args.index('-dt') + 1]))
		elif '-h' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-p <sequence of projection files> - this param must be the first\n -pr <range of files: "proj1.xvg-proj100.xvg">\n-o <timescales file> (*.xvg)\n\
-ln <number of lipids>\n -dt <timestep in (ns)>')


	elif args[0] == 'autot':
		main = autocorr.main
		if '-p' in args and '-ln' in args and '-dt' in args and '-o' in args and '-pr' not in args:
			filenames = [args[i] for i in range(args.index('-p') + 1, args.index('-o'))]
			main(filenames, int(args[args.index('-ln') + 1]), float(args[args.index('-dt') + 1]), args[args.index('-o') + 1])
		elif '-p' in args and '-ln' in args and '-dt' in args and '-o' not in args and '-pr' not in args:
			# print('No output file supplied. Data will be written in "autocorr_relaxtime_vs_PC.xvg"')
			filenames = [args[i] for i in range(args.index('-p') + 1, min(args.index('-dt'), args.index('-ln')))]
			main(filenames, int(args[args.index('-ln') + 1]), float(args[args.index('-dt') + 1]))
		elif '-pr' in args and '-ln' in args and '-dt' in args and '-o' in args and '-p' not in args:
			files = args[args.index('-pr') + 1]
			file_start = files[:files.find('-')]
			file_end = files[files.find('-') + 1:]
			start = int(file_start[file_start.rfind('_') + 1:file_start.rfind('.')])
			end = int(file_end[file_end.rfind('_') + 1:file_end.rfind('.')])
			file_mask = file_start[:file_start.rfind('_')] + file_start[file_start.rfind('.'):]
			filenames = [file_mask[:file_mask.find('.')] + "_" + str(i) + file_mask[file_mask.find('.'):] for i in range(start, end + 1)]
			main(filenames, int(args[args.index('-ln') + 1]), float(args[args.index('-dt') + 1]), args[args.index('-o') + 1])
		elif '-pr' in args and '-ln' in args and '-dt' in args and '-o' not in args and '-p' not in args: 
			files = args[args.index('-pr') + 1]
			file_start = files[:files.find('-')]
			file_end = files[files.find('-') + 1:]
			start = int(file_start[file_start.rfind('_') + 1:file_start.rfind('.')])
			end = int(file_end[file_end.rfind('_') + 1:file_end.rfind('.')])
			file_mask = file_start[:file_start.rfind('_')] + file_start[file_start.rfind('.'):]
			filenames = [file_mask[:file_mask.find('.')] + "_" + str(i) + file_mask[file_mask.find('.'):] for i in range(start, end + 1)]
			main(filenames, int(args[args.index('-ln') + 1]), float(args[args.index('-dt') + 1]))
		elif '-h' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-p <sequence of projection files> - this param must be the first\n -pr <range of files: "proj1.xvg-proj100.xvg">\n-o <timescales file> (*.xvg)\n\
-ln <number of lipids>\n-dt <timestep in (ns)>')

	elif args[0] == 'eigenvecdot':
		main = eigenvecdot.main
		if '-evec' in args:
			main(args[args.index('-evec') + 1], args[args.index('-evec') + 2])
		elif '-h' not in args and '-help' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-evec <two files with eigenvectors (example: -evec eigenvec1.xvg eigenvec2.xvg)>')


	elif args[0] == 'conspace':
		main = lipic.main
		if '-stride' in args:
			stride = int(args[args.index('-stride') + 1])
		else:
			stride = 10000
		if '-om' in args:
			mot_out = args[args.index('-om') + 1]
		else:
			mot_out = 'conform.pdb'
		if '-f' in args and '-t' in args:
			main(args[args.index('-f') + 1], args[args.index('-t') + 1], stride, mot_out)
		elif '-h' not in args and '-help' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-f <trajectory file> (file format *.xtc)\n-t <topology file> (any file with topology)\n -stride <<positive integer; step of reading frames>\n -om <output file with conformations>')


	elif args[0] == 'pearson':
		main = pearson.main
		if '-cov1' in args and '-cov2' in args:
			main(args[args.index('-cov1') + 1], args[args.index('-cov2') + 1])
		elif '-h' not in args and '-help' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-cov1, -cov2 - 2 files with covariance matrices')

	elif args[0] == 'splitproj':
		main = proj_parse.main
		if '-p' in args:
			main(args[args.index('-p') + 1])
		elif '-h' not in args and '-help' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-p - projection file')

	elif args[0] == 'motion':
		main = visualizing.main
		if '-p' in args and '-npc' in args and '-aver' in args and '-ievec' in args:
			main(args[args.index('-p') + 1], args[args.index('-npc') + 1], args[args.index('-aver') + 1], args[args.index('-ievec') + 1])
		else:
			print('-p - projection file\n -npc - number of principal component\n -aver - average structure\n -ievec - file with eigenvectors')

	elif args[0] == 'projdistm':
		main = projdistm.main
		if '-files' in args and len(args) > 2:
			main(args[args.index('-files') + 1:])
		elif '-h' not in args and '-help' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-files <sequence of files with projection>')


	elif args[0] == 'reltime':
		main = reltime.main
		if '-eval' in args and '-time1' in args and '-time2' in args:
			main(args[args.index('-eval') + 1],args[args.index('-time1') + 1],args[args.index('-time2') + 1])
		elif '-h' not in args and '-help' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-eval - file with eigenvalues \n-time1 and -time2 - 2 files related to different trajectories that contains relaxation time for autocorrelations')


	elif args[0] == 'combtrajs':
		main = combtraj.main
		if '-fs' in args:
			i = args.index('-fs')
			length = len(args[(i+1):])
			omaewa = []
			for j in range(i+1, len(args), 2):
				omaewa.append((args[j], args[j + 1]))
			main(omaewa)
		elif '-h' not in args and '-help' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-fs - input files with trajectories and topologies in the right order: \n-fs concatenated_1.xtc average_1.pdb concatenated_2.xtc average_2.pdb')


	elif args[0] == 'timescalespic':
		main = timescalespic.main
		if '-file1' in args and '-file2' in args and '-type' in args and '-time' in args:
			main(args[args.index('-file1') + 1],args[args.index('-file2') + 1],args[args.index('-type') + 1],args[args.index('-time') + 1])
		elif '-h' not in args and '-help' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-file1 and -file2 - input files with timescales\n-type "kss" or "auto" for kss or autocorrelation data\n-time "t1" or "t2" for decreasing in e or e^2 times)')


	elif args[0] == 'eigenvals':
		main = eigenvalues.main
		if '-cumulative' in args:
			cumul = args[args.index('-cumulative') + 1]
		if cumul == "True" or cumul == "1":
			cumul = True
		else:
			cumul = False
		if '-ieval' in args:
			main(args[args.index('-ieval') + 1], cumul)
		elif '-h' not in args and '-help' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-ieval - input file with eigevalues (file format *.xvg)\n-cumulative - cumulative sum of eigenvalues (optional <False>)')


	elif args[0] == 'help' or args[0] == '-h' or args[0] == '-help':
		print("'concat' - create concatenated trajectory\n\
'covar' - principal component analysis\n\
'project' - calculating projections\n\
'projdist' - probability density of single trajectory projection\n\
'projdistm' - probability density of two trajectories projections\n\
'ksst' - Kolmogorov-Smirnov convergence\n\
'autot' - Autocorrelation decay\n\
'eigenvecdot' - scalar product of eigenvectors from different trajectories\n\
'conspace' - conformational space of lipids in trajectory\n\
'pearson' - Pearson coefficient for covariance matrices from different trajectories\n\
'splitproj' - split files with all projections into files with projections for single PC\n\
'motion' - demonstrates the motion along PC.\n\
'eigenvals' - picture of eigenvalues of covariance matrix or their cumulative sum\n\
'combtrajs' - combine 2 concatenated trajectories into one associated trajectory\n\
'reltime' - comparison of the characteristic timescales for KSS or autocorrelation\n\
'timescalespic' - picture for joint analysis of the timescales of two trajectories\n\
 Use any of this options.")

	else:
		print('Use -h or help for more information.')


if __name__ == '__main__':
	args = sys.argv[1:]
	if len(args) == 0:
		print('Use -h or help for more information.')
	else:
		main(args)
