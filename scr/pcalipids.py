import PCA
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
			print(' -f <trajectory file> (file format *.xtc, *trr)\n -t <topology file> (any file with topology)\n -r <reference traj file> (any topology file). If not supplied, \
the first frame of trajectory will be used for alignment\n -l <lipid type> (example: -l DPPC)\n -stride <positive integer; step of reading frames>\n \
-sf <time in ps; number to determine from which frame to read the trajectory>\n -oc <output trajectory file>\n -oa <output topology file>\n')

	
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
		if '-op' in args:
			proj_file = args[args.index('-op') + 1]
		else:
			proj_file = None
		if '-f' in args and '-t' in args and '-first' in args and '-last' in args:
			main(args[args.index('-f') + 1], args[args.index('-t') + 1], first_PC = args[args.index('-first') + 1], last_PC = args[args.index('-last') + 1], val_file = val_file, vec_file = vec_file, cov_file = cov_file, proj_file = proj_file)
		elif '-f' in args and '-t' in args and '-first' not in args and '-last' not in args:
			main(args[args.index('-f') + 1], args[args.index('-t') + 1], val_file = val_file, vec_file = vec_file, cov_file = cov_file, proj_file = proj_file, first_PC = None, last_PC = None)
		elif '-f' in args and '-t' in args and '-first' in args and '-last' not in args:
			main(args[args.index('-f') + 1], args[args.index('-t') + 1], first_PC = args[args.index('-first') + 1], val_file = val_file, vec_file = vec_file, cov_file = cov_file, proj_file = proj_file, last_PC = None)
		elif '-f' in args and '-t' in args and '-first' not in args and '-last' in args:
			main(args[args.index('-f') + 1], args[args.index('-t') + 1], last_PC = args[args.index('-last') + 1], val_file = val_file, vec_file = vec_file, cov_file = cov_file, proj_file = proj_file, first_PC =None)
		elif '-h' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-f <trajectory file> (file format *.xtc)\n-t <topology file> (any file with topology)\n -first <first PC> -last <last PC> \n -oeval <output file with eigenvalues>\n\
				-oevec <output file with eigenvectors>\n -ocov <output file with covariance matrix>\n -op <output file with projections>.')

	
	elif args[0] == 'projdist':
		main = proj_dist_s.main
		if '-p' in args and '-first' in args and '-last' in args:
			main(args[args.index('-p') + 1], args[args.index('-first') + 1], args[args.index('-last') + 1])
		elif '-p' in args and '-first' not in args and '-last' not in args:
			main(args[args.index('-p') + 1], 1, 3)
		elif '-h' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-p <projection file> (file format *.xvg)\n -first <first projection> -last <last projection> (int format). \nIf not supplied, the first 3 projections will be analyzed.')

	
	elif args[0] == 'ksst':
		main = kss_time_scale.main
		if '-p' in args and '-ln' in args and '-o' in args and '-pr' not in args:
			start = args.index('-p')
			end = min(args.index('-o'), args.index('-ln'))
			filenames = [args[i] for i in range(start + 1, end)]
			main(filenames, int(args[args.index('-ln') + 1]), args[args.index('-o') + 1])
		elif '-p' in args and '-ln' in args and '-o' not in args and '-pr' not in args:
			print('No output file supplied. Data will be written in "KSS_relaxtime_vs_PC.xvg"')
			filenames = [args[i] for i in range(args.index('-p') + 1, args.index('-o'))]
			main(filenames, int(args[args.index('-ln') + 1]))
		elif '-pr' in args and '-ln' in args and '-o' in args and '-p' not in args:
			files = args[args.index('-pr') + 1]
			file_start = files[:files.find('-')]
			file_end = files[files.find('-') + 1:]
			start = int(''.join(filter(lambda x: x.isdigit(), file_start)))
			end = int(''.join(filter(lambda x: x.isdigit(), file_end)))
			file_mask = ''.join(filter(lambda x: not x.isdigit(), file_start))
			filenames = [file_mask[:file_mask.find('.')] + str(i) + file_mask[file_mask.find('.'):] for i in range(start, end + 1)]
			main(filenames, int(args[args.index('-ln') + 1]), args[args.index('-o') + 1])
		elif '-pr' in args and '-ln' in args and '-o' not in args and '-p' not in args: 
			files = args[args.index('-pr') + 1]
			file_start = files[:files.find('-')]
			file_end = files[files.find('-') + 1:]
			start = int(''.join(filter(lambda x: x.isdigit(), file_start)))
			end = int(''.join(filter(lambda x: x.isdigit(), file_end)))
			file_mask = ''.join(filter(lambda x: not x.isdigit(), file_start))
			filenames = [file_mask[:file_mask.find('.')] + str(i) + file_mask[file_mask.find('.'):] for i in range(start, end + 1)]
			main(filenames, int(args[args.index('-ln') + 1]))
		elif '-h' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-p <sequence of projection files> - this param must be the first\n -pr <range of files: "proj1.xvg-proj100.xvg">\n-t <topology file> (any file with topology)\n -ln <number of lipids>\n-r <reference traj file> (any file with 1 ref frame). If not supplied, \
				the first frame of trajectory will be used for alignment.')


	elif args[0] == 'autot':
		main = autocorr.main
		if '-p' in args and '-o' in args and '-pr' not in args:
			filenames = [args[i] for i in range(args.index('-p') + 1, args.index('-o'))]
			main(filenames, args[args.index('-o') + 1])
		elif '-p' in args and '-o' not in args and '-pr' not in args:
			print('No output file supplied. Data will be written in "autocorr_relaxtime_vs_PC.xvg"')
			filenames = [args[i] for i in range(args.index('-p') + 1, args.index('-o'))]
			main(filenames)
		elif '-pr' in args and '-o' in args and '-p' not in args:
			files = args[args.index('-pr') + 1]
			file_start = files[:files.find('-')]
			file_end = files[files.find('-') + 1:]
			start = int(''.join(filter(lambda x: x.isdigit(), file_start)))
			end = int(''.join(filter(lambda x: x.isdigit(), file_end)))
			file_mask = ''.join(filter(lambda x: not x.isdigit(), file_start))
			filenames = [file_mask[:file_mask.find('.')] + str(i) + file_mask[file_mask.find('.'):] for i in range(start, end + 1)]
			main(filenames, args[args.index('-o') + 1])
		elif '-pr' in args and '-o' not in args and '-p' not in args: 
			files = args[args.index('-pr') + 1]
			file_start = files[:files.find('-')]
			file_end = files[files.find('-') + 1:]
			start = int(''.join(filter(lambda x: x.isdigit(), file_start)))
			end = int(''.join(filter(lambda x: x.isdigit(), file_end)))
			file_mask = ''.join(filter(lambda x: not x.isdigit(), file_start))
			filenames = [file_mask[:file_mask.find('.')] + str(i) + file_mask[file_mask.find('.'):] for i in range(start, end + 1)]
			main(filenames)
		elif '-h' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-p <sequence of projection files>\n -pr <range of files: "proj1.xvg-proj100.xvg">\n-t <topology file> (any file with topology)\n-r <reference traj file> (any file with 1 ref frame). If not supplied, \
				the first frame of trajectory will be used for alignment.')

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


	elif args[0] == 'projdist':
		main = proj_dist_s.main
		if '-p' in args and '-first' in args and '-last' in args:
	 		main(args[args.index('-p') + 1], args[args.index('-first') + 1], args[args.index('-last') + 1])
		elif '-p' in args and '-first' not in args and '-last' not in args:
	 		main(args[args.index('-p') + 1], 1, 3)
		elif '-h' not in args:
	 		print('Missing parameters, try -h for flags\n')
		else:
			print('-p <projection file> (file format *.xvg)\n -fisrt <first projection> -last <last projection> (int format). \nIf not supplied, the first 3 projections will be analyzed.')

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
		if '-p' in args and '-npc' in args and '-aver' in args and '-e' in args:
			main(args[args.index('-p') + 1], args[args.index('-npc') + 1], args[args.index('-aver') + 1], args[args.index('-e') + 1])
		else:
			print('-p - projection file\n -npc - number of principal component\n -aver - average structure\n -e - file with eigenvectors')


	elif args[0] == 'help' or args[0] == '-h' or args[0] == '-help':
		print("'concat' - create concatenated trajectory\n\
'covar' - principal component analysis\n\
'ksst' - Kolmogorov-Smirnov convergence\n\
'autot' - Autocorrelation decay\n\
'eigenvecdot' - scalar product of eigenvectors from different trajectories\
'conspace' - conformational space of lipids in trajectory\
'pearson' - Pearson coefficient for covariance matrices from different trajectories\
'splitproj - split files with all projections into files with projections for single PC\
'motion' - demonstrates the motion along PC.\n\
Use any of this options.")

	else:
		print('Use -h or help for more information.')


if __name__ == '__main__':
	args = sys.argv[1:]
	if len(args) == 0:
		print('Use -h or help for more information.')
	else:
		main(args)