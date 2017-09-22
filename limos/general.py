import PCA
import con_traj
import parse_proj_file
import sys


def main(args):
	if  args[0] == 'concat':
		if '-f' in args and '-t' in args and '-r' in args:
			con_traj.main(args[args.index('-f') + 1], args[args.index('-t') + 1], args[args.index('-r') + 1])
		elif '-f' in args and '-t' in args and '-r' not in args:
			print('No reference file supplied. The first frame of trajectory will be used for alignment.')
			con_traj.main(args[args.index('-f') + 1], args[args.index('-t') + 1])
		elif '-h' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-f <trajectory file> (file format *.xtc)\n-t <topology file> (any file with topology)\n-r <reference traj file> (any file with 1 ref frame). If not supplied, \
				the first frame of trajectory will be used for alignment.')

	elif args[0] == 'covar':
		if '-f' in args and '-t' in args and '-first' in args and '-last' in args:
			PCA.main(args[args.index('-f') + 1], args[args.index('-t') + 1], args[args.index('-first') + 1], args[args.index('-last') + 1])
		elif '-f' in args and '-t' in args and '-first' not in args and '-last' not in args:
			PCA.main(args[args.index('-f') + 1], args[args.index('-t') + 1], None, None)
		elif '-f' in args and '-t' in args and '-first' in args and '-last' not in args:
			PCA.main(args[args.index('-f') + 1], args[args.index('-t') + 1], args[args.index('-first') + 1], None)
		elif '-f' in args and '-t' in args and '-first' not in args and '-last' in args:
			PCA.main(args[args.index('-f') + 1], args[args.index('-t') + 1], None, args[args.index('-last') + 1])
		elif '-h' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-f <trajectory file> (file format *.xtc)\n-t <topology file> (any file with topology)\n -first <first PC> -last <last PC>.')

	elif args[0] == 'anaproj':
		if '-p' in args and '-first' in args and '-last' in args:
			parse_proj_file.main(args[args.index('-p') + 1], args[args.index('-first') + 1], args[args.index('-last') + 1])
		elif '-p' in args and '-first' not in args and '-last' not in args:
			parse_proj_file.main(args[args.index('-p') + 1], 1, 3)
		elif '-h' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-p <projection file> (file format *.xvg)\n -fisrt <first projection> -last <last projection> (int format). \nIf not supplied, the first 3 projections will be analyzed.')

	elif args[0] == 'all':
		if '-f' in args and '-t' in args and '-r' in args:
			traj_, top_ = con_traj.main(args[args.index('-f') + 1], args[args.index('-t') + 1], args[args.index('-r') + 1])
		elif '-f' in args and '-t' in args and '-r' not in args:
			print('No reference file supplied. The first frame of trajectory will be used for alignment.')
			traj_, top_ = con_traj.main(args[args.index('-f') + 1], args[args.index('-t') + 1])
		elif '-h' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-f <trajectory file> (file format *.xtc)\n-t <topology file> (any file with topology)\n-r <reference traj file> (any file with 1 ref frame). If not supplied, \
				the first frame of trajectory will be used for alignment.')

		if '-first' in args and '-last' in args:
			proj_, first_, last_ = PCA.main(traj_, top_, args[args.index('-first') + 1], args[args.index('-last') + 1])
		elif '-first' not in args and '-last' not in args:
			proj_, first_, last_ = PCA.main(traj_, top_, args[args.index('-f') + 1], args[args.index('-t') + 1], None, None)
		elif '-first' in args and '-last' not in args:
			proj_, first_, last_ = PCA.main(traj_, top_, args[args.index('-f') + 1], args[args.index('-t') + 1], args[args.index('-first') + 1], None)
		elif '-first' not in args and '-last' in args:
			proj_, first_, last_ = PCA.main(traj_, top_, args[args.index('-f') + 1], args[args.index('-t') + 1], None, args[args.index('-last') + 1])
		elif '-h' not in args:
			print('Missing parameters, try -h for flags\n')
		else:
			print('-f <trajectory file> (file format *.xtc)\n-t <topology file> (any file with topology)\n -first <first PC> -last <last PC>.')

		parse_proj_file.main(proj_, first_, last_)

	elif args[0] == 'help' or args[0] == '-h' or args[0] == '-help':
		print("'concat' - create concatenated trajectory\n\
'covar' - principal component analysis\n\
'anaproj' - projections distribution\n\
'all' - all steps\n\
Use any of this options.")

	else:
		print('Use -h or help for more information.')


if __name__ == '__main__':
	args = sys.argv[1:]
	if len(args) == 0:
		print('Use -h or help for more information.')
	else:
		main(args)