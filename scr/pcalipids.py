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
import argparse

# Basic option class
# Needed to select correct function to continue
class Option:
    def __init__(self,func=str,fname=str,optList=list,description=""):
        self.func        = func
        self.fname       = fname
        self.optList     = optList
        self.description = description

# Function option class
# Needed to parse the input parameters for function
class OptionF:
    def __init__(self,func=str,num=1,default=None,description=""):
        self.func        = func
        self.num         = num
        self.value       = default
        self.description = description
    def __nonzero__(self): 
        return self.value != None
    def __str__(self):
        return self.value and str(self.value) or ""
    def setvalue(self,v):
        if len(v) == 1:
            self.value = self.func(v[0])
        else:
            self.value = [ self.func(i) for i in v ]
    
def main(args):

	optionsConcat = [
	# options for concat feature
	"Input/output options for concat feature",
	("-f", OptionF(str, 1, None, "Input XTC or TRR file")),
	("-t", OptionF(str, 1, None, "Input topology PDB of GRO file")),
	("-stride", OptionF(int, 1, 1, "Only read every Nr-th frame")),
	("-sf", OptionF(int, 1, 0, "First frame (ps) to read from trajectory")),
	("-ef", OptionF(int, 1, -1, "Last frame (ps) to read from trajectory")),
	("-oc", OptionF(str, 1, "concatenated.xtc", "Output concatenated trajectory file")),
	("-oa", OptionF(str, 1, "average.pdb", "Output average structure")),
	("-r", OptionF(str, 1, None, "Input reference PDB of GRO file. If not supplied, \
		the structure for the first lipid and the first frame is used for alignment")),
	("-l", OptionF(str, 1, None, "Lipid type"))		
	]

	optionsCovar = [
	# options for covar feature
	"Input/output options for covar feature",
	("-f", OptionF(str, 1, None, "Input XTC or TRR concatenated trajectory")),
	("-t", OptionF(str, 1, None, "Input topology PDB of GRO file")),
	("-oeval", OptionF(str, 1, "eigenval.xvg", "Eigenvalue output file")),
	("-oevec", OptionF(str, 1, "eigenvec.xvg", "Eigenvector output file")),
	("-ocov", OptionF(str, 1, "covar.dat", "Covariance matrix output file")),
	("-invertPC1", OptionF(int, 1, 0, "Invert the distribution for the first PC or not\n\
		[default] 0 - do not invert\n\
		1 - invert"))
	]

	optionsProject = [
	# options for project feature
	"Input/output options for project feature",
	("-f", OptionF(str, 1, None, "Input XTC or TRR concatenated trajectory")),
	("-t", OptionF(str, 1, None, "Input topology PDB of GRO file")),
	("-ia", OptionF(str, 1, None, "PDB of GRO file for average structure")),
	("-ievec", OptionF(str, 1, None, "Eigenvector input file")),
	("-first", OptionF(int, 1, 1, "First PC for projection (default: 1)")),
	("-last", OptionF(int, 1, 10, "First PC for projection (default: 10)")),
	("-op", OptionF(str, 1, "proj.xvg", "Output projection files name\n\
		Do not use '-' or '.' symbols in the projection file names"))
	]

	optionsProjDist = [
	# options for projdist feature
	"Input/output options for projdist feature",
	("-p", OptionF(str, 1, None, "Input projection file")),
	("-pr", OptionF(str, 1, None, "Range of input projection files: \n\
	example: -pr proj_1.xvg-proj_10.xvg"))
	]

	options = [
	# list of all available features
	"List of procedures",
	("concat", Option(str,"con_traj",optionsConcat,\
		"Concatenate trajectories of individual lipids")),
	("covar", Option(str,"PCA",optionsCovar,\
		"Perform PCA on aligned concatenated lipid trajectory")),
	("project", Option(str,"project",optionsProject,\
		"Project concatenated lipid trajectory on the calculated PCs")),
	("projdist", Option(str,"proj_dist_s",optionsProjDist,\
		"Plot projection distributions for selected PCs"))
	]

	desc = "\nPCAlipids is a software for analysis of lipid molecule conformations\n"

	# If the user asks for help: pcalipids.py -h
	if args[0] == '-h' or args[0] == '--help':
		print("\n",__file__) # print executable file name
		print(desc) # print description !!! we need to create a good one
		for thing in options: # print all options
			print(type(thing) != str and "%10s: %s"%(thing[0],thing[1].description) or thing)
		print()
		sys.exit()

	# Create a dictionary from options
	# All strings are ignored (they descriptive)
	options = dict([i for i in options if not type(i) == str])

	# Redirect main to the correct import class
	main = globals()[options[args[0]].fname].main
	optionsF = options[args[0]].optList

	# We don't need 0 argument of args any more
	args = args[1:]

	# if the user asks for help: pcalipids.py concat -h
	if args[0] == '-h' or args[0] == '--help':
		print("\n",__file__) # print executable file name
		print(desc) # print description !!! we need to create a good one
		for thing in optionsF: # print all options for selected feature
			print(type(thing) != str and "%10s: %s"%(thing[0],thing[1].description) or thing)
		print()
		sys.exit()

	# Create a dictionary with possible arguments for the selected feature
	optionsF = dict([i for i in optionsF if not type(i) == str])

	# Parsing argumens and setting values
	while args:
		ar = args.pop(0) # choose argument
		optionsF[ar].setvalue([args.pop(0) for i in range(optionsF[ar].num)]) # set value

	# Pass all the parameters to function
	params=list(v.value for v in optionsF.values())
	main(*params)

# 	elif args[0] == 'ksst':
# 		main = kss_time_scale.main
# 		if '-p' in args and '-ln' in args and '-dt' in args and '-o' in args and '-pr' not in args:
# 			start = args.index('-p')
# 			end = min(args.index('-o'), args.index('-ln'))
# 			filenames = [args[i] for i in range(start + 1, end)]
# 			main(filenames, int(args[args.index('-ln') + 1]), float(args[args.index('-dt') + 1]), args[args.index('-o') + 1])
# 		elif '-p' in args and '-ln' in args and '-dt' in args and '-o' not in args and '-pr' not in args:
# 			# print('No output file supplied. Data will be written in "KSS_relaxtime_vs_PC.xvg"')
# 			filenames = [args[i] for i in range(args.index('-p') + 1, min(args.index('-dt'), args.index('-ln')))]
# 			main(filenames, int(args[args.index('-ln') + 1]), float(args[args.index('-dt') + 1]))
# 		elif '-pr' in args and '-ln' in args and '-dt' in args and '-o' in args and '-p' not in args:
# 			files = args[args.index('-pr') + 1]
# 			file_start = files[:files.find('-')]
# 			file_end = files[files.find('-') + 1:]
# 			start = int(file_start[file_start.rfind('_') + 1:file_start.rfind('.')])
# 			end = int(file_end[file_end.rfind('_') + 1:file_end.rfind('.')])
# 			file_mask = file_start[:file_start.rfind('_')] + file_start[file_start.rfind('.'):]
# 			filenames = [file_mask[:file_mask.find('.')] + "_" + str(i) + file_mask[file_mask.find('.'):] for i in range(start, end + 1)]
# 			main(filenames, int(args[args.index('-ln') + 1]), float(args[args.index('-dt') + 1]), args[args.index('-o') + 1])
# 		elif '-pr' in args and '-ln' in args and '-dt' in args and '-o' not in args and '-p' not in args: 
# 			files = args[args.index('-pr') + 1]
# 			file_start = files[:files.find('-')]
# 			file_end = files[files.find('-') + 1:]
# 			start = int(file_start[file_start.rfind('_') + 1:file_start.rfind('.')])
# 			end = int(file_end[file_end.rfind('_') + 1:file_end.rfind('.')])
# 			file_mask = file_start[:file_start.rfind('_')] + file_start[file_start.rfind('.'):]
# 			filenames = [file_mask[:file_mask.find('.')] + "_" + str(i) + file_mask[file_mask.find('.'):] for i in range(start, end + 1)]
# 			main(filenames, int(args[args.index('-ln') + 1]), float(args[args.index('-dt') + 1]))
# 		elif '-h' not in args:
# 			print('Missing parameters, try -h for flags\n')
# 		else:
# 			print('-p <sequence of projection files> - this param must be the first\n -pr <range of files: "proj1.xvg-proj100.xvg">\n-o <timescales file> (*.xvg)\n\
# -ln <number of lipids>\n -dt <timestep in (ns)>')


# 	elif args[0] == 'autot':
# 		main = autocorr.main
# 		if '-p' in args and '-ln' in args and '-dt' in args and '-o' in args and '-pr' not in args:
# 			filenames = [args[i] for i in range(args.index('-p') + 1, args.index('-o'))]
# 			main(filenames, int(args[args.index('-ln') + 1]), float(args[args.index('-dt') + 1]), args[args.index('-o') + 1])
# 		elif '-p' in args and '-ln' in args and '-dt' in args and '-o' not in args and '-pr' not in args:
# 			# print('No output file supplied. Data will be written in "autocorr_relaxtime_vs_PC.xvg"')
# 			filenames = [args[i] for i in range(args.index('-p') + 1, min(args.index('-dt'), args.index('-ln')))]
# 			main(filenames, int(args[args.index('-ln') + 1]), float(args[args.index('-dt') + 1]))
# 		elif '-pr' in args and '-ln' in args and '-dt' in args and '-o' in args and '-p' not in args:
# 			files = args[args.index('-pr') + 1]
# 			file_start = files[:files.find('-')]
# 			file_end = files[files.find('-') + 1:]
# 			start = int(file_start[file_start.rfind('_') + 1:file_start.rfind('.')])
# 			end = int(file_end[file_end.rfind('_') + 1:file_end.rfind('.')])
# 			file_mask = file_start[:file_start.rfind('_')] + file_start[file_start.rfind('.'):]
# 			filenames = [file_mask[:file_mask.find('.')] + "_" + str(i) + file_mask[file_mask.find('.'):] for i in range(start, end + 1)]
# 			main(filenames, int(args[args.index('-ln') + 1]), float(args[args.index('-dt') + 1]), args[args.index('-o') + 1])
# 		elif '-pr' in args and '-ln' in args and '-dt' in args and '-o' not in args and '-p' not in args: 
# 			files = args[args.index('-pr') + 1]
# 			file_start = files[:files.find('-')]
# 			file_end = files[files.find('-') + 1:]
# 			start = int(file_start[file_start.rfind('_') + 1:file_start.rfind('.')])
# 			end = int(file_end[file_end.rfind('_') + 1:file_end.rfind('.')])
# 			file_mask = file_start[:file_start.rfind('_')] + file_start[file_start.rfind('.'):]
# 			filenames = [file_mask[:file_mask.find('.')] + "_" + str(i) + file_mask[file_mask.find('.'):] for i in range(start, end + 1)]
# 			main(filenames, int(args[args.index('-ln') + 1]), float(args[args.index('-dt') + 1]))
# 		elif '-h' not in args:
# 			print('Missing parameters, try -h for flags\n')
# 		else:
# 			print('-p <sequence of projection files> - this param must be the first\n -pr <range of files: "proj1.xvg-proj100.xvg">\n-o <timescales file> (*.xvg)\n\
# -ln <number of lipids>\n-dt <timestep in (ns)>')

# 	elif args[0] == 'eigenvecdot':
# 		main = eigenvecdot.main
# 		if '-evec' in args:
# 			main(args[args.index('-evec') + 1], args[args.index('-evec') + 2])
# 		elif '-h' not in args and '-help' not in args:
# 			print('Missing parameters, try -h for flags\n')
# 		else:
# 			print('-evec <two files with eigenvectors (example: -evec eigenvec1.xvg eigenvec2.xvg)>')


# 	elif args[0] == 'conspace':
# 		main = lipic.main
# 		if '-stride' in args:
# 			stride = int(args[args.index('-stride') + 1])
# 		else:
# 			stride = 10000
# 		if '-om' in args:
# 			mot_out = args[args.index('-om') + 1]
# 		else:
# 			mot_out = 'conform.pdb'
# 		if '-f' in args and '-t' in args:
# 			main(args[args.index('-f') + 1], args[args.index('-t') + 1], stride, mot_out)
# 		elif '-h' not in args and '-help' not in args:
# 			print('Missing parameters, try -h for flags\n')
# 		else:
# 			print('-f <trajectory file> (file format *.xtc)\n-t <topology file> (any file with topology)\n -stride <<positive integer; step of reading frames>\n -om <output file with conformations>')


# 	elif args[0] == 'pearson':
# 		main = pearson.main
# 		if '-cov1' in args and '-cov2' in args:
# 			main(args[args.index('-cov1') + 1], args[args.index('-cov2') + 1])
# 		elif '-h' not in args and '-help' not in args:
# 			print('Missing parameters, try -h for flags\n')
# 		else:
# 			print('-cov1, -cov2 - 2 files with covariance matrices')

# 	elif args[0] == 'splitproj':
# 		main = proj_parse.main
# 		if '-p' in args:
# 			main(args[args.index('-p') + 1])
# 		elif '-h' not in args and '-help' not in args:
# 			print('Missing parameters, try -h for flags\n')
# 		else:
# 			print('-p - projection file')

# 	elif args[0] == 'motion':
# 		main = visualizing.main
# 		if '-p' in args and '-npc' in args and '-aver' in args and '-ievec' in args:
# 			main(args[args.index('-p') + 1], args[args.index('-npc') + 1], args[args.index('-aver') + 1], args[args.index('-ievec') + 1])
# 		else:
# 			print('-p - projection file\n -npc - number of principal component\n -aver - average structure\n -ievec - file with eigenvectors')

# 	elif args[0] == 'projdistm':
# 		main = projdistm.main
# 		if '-files' in args and len(args) > 2:
# 			main(args[args.index('-files') + 1:])
# 		elif '-h' not in args and '-help' not in args:
# 			print('Missing parameters, try -h for flags\n')
# 		else:
# 			print('-files <sequence of files with projection>')


# 	elif args[0] == 'reltime':
# 		main = reltime.main
# 		if '-eval' in args and '-time1' in args and '-time2' in args:
# 			main(args[args.index('-eval') + 1],args[args.index('-time1') + 1],args[args.index('-time2') + 1])
# 		elif '-h' not in args and '-help' not in args:
# 			print('Missing parameters, try -h for flags\n')
# 		else:
# 			print('-eval - file with eigenvalues \n-time1 and -time2 - 2 files related to different trajectories that contains relaxation time for autocorrelations')


# 	elif args[0] == 'combtrajs':
# 		main = combtraj.main
# 		if '-fs' in args:
# 			i = args.index('-fs')
# 			length = len(args[(i+1):])
# 			omaewa = []
# 			for j in range(i+1, len(args), 2):
# 				omaewa.append((args[j], args[j + 1]))
# 			main(omaewa)
# 		elif '-h' not in args and '-help' not in args:
# 			print('Missing parameters, try -h for flags\n')
# 		else:
# 			print('-fs - input files with trajectories and topologies in the right order: \n-fs concatenated_1.xtc average_1.pdb concatenated_2.xtc average_2.pdb')


# 	elif args[0] == 'timescalespic':
# 		main = timescalespic.main
# 		if '-file1' in args and '-file2' in args and '-type' in args and '-time' in args:
# 			main(args[args.index('-file1') + 1],args[args.index('-file2') + 1],args[args.index('-type') + 1],args[args.index('-time') + 1])
# 		elif '-h' not in args and '-help' not in args:
# 			print('Missing parameters, try -h for flags\n')
# 		else:
# 			print('-file1 and -file2 - input files with timescales\n-type "kss" or "auto" for kss or autocorrelation data\n-time "t1" or "t2" for decreasing in e or e^2 times)')


# 	elif args[0] == 'eigenvals':
# 		main = eigenvalues.main
# 		if '-cumulative' in args:
# 			cumul = args[args.index('-cumulative') + 1]
# 		if cumul == "True" or cumul == "1":
# 			cumul = True
# 		else:
# 			cumul = False
# 		if '-ieval' in args:
# 			main(args[args.index('-ieval') + 1], cumul)
# 		elif '-h' not in args and '-help' not in args:
# 			print('Missing parameters, try -h for flags\n')
# 		else:
# 			print('-ieval - input file with eigevalues (file format *.xvg)\n-cumulative - cumulative sum of eigenvalues (optional <False>)')


# 	elif args[0] == 'help' or args[0] == '-h' or args[0] == '-help':
# 		print("'concat' - create concatenated trajectory\n\
# 'covar' - principal component analysis\n\
# 'project' - calculating projections\n\
# 'projdist' - probability density of single trajectory projection\n\
# 'projdistm' - probability density of two trajectories projections\n\
# 'ksst' - Kolmogorov-Smirnov convergence\n\
# 'autot' - Autocorrelation decay\n\
# 'eigenvecdot' - scalar product of eigenvectors from different trajectories\n\
# 'conspace' - conformational space of lipids in trajectory\n\
# 'pearson' - Pearson coefficient for covariance matrices from different trajectories\n\
# 'splitproj' - split files with all projections into files with projections for single PC\n\
# 'motion' - demonstrates the motion along PC.\n\
# 'eigenvals' - picture of eigenvalues of covariance matrix or their cumulative sum\n\
# 'combtrajs' - combine 2 concatenated trajectories into one associated trajectory\n\
# 'reltime' - comparison of the characteristic timescales for KSS or autocorrelation\n\
# 'timescalespic' - picture for joint analysis of the timescales of two trajectories\n\
#  Use any of this options.")

# 	else:
# 		print('Use -h or help for more information.')
###

if __name__ == '__main__':
	args = sys.argv[1:]
	main(args)
