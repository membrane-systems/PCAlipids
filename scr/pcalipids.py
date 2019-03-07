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
	("-f", OptionF(str, 1, None, "Input trajectory file (.xtc, .trr, ...)")),
	("-t", OptionF(str, 1, None, "Input topology file (.pdb, .gro, ...)")),
	("-stride", OptionF(int, 1, 1, "Only read every Nth frame")),
	("-sf", OptionF(int, 1, 0, "First frame (ps) to read from trajectory")),
	("-ef", OptionF(int, 1, -1, "Last frame (ps) to read from trajectory")),
	("-oc", OptionF(str, 1, "concatenated.xtc", "Output concatenated trajectory file")),
	("-oa", OptionF(str, 1, "average.pdb", "Output average structure")),
	("-r", OptionF(str, 1, None, "Input reference file (.pdb, .gro). If not supplied, \
the structure of the first lipid and the first frame is used for alignment")),
	("-l", OptionF(str, 1, None, "Lipid type"))		
	]

	optionsConspace = [
	# options for conspace feature
	"Input/output options for conspace feature",
	("-f", OptionF(str, 1, None, "Input XTC or TRR file")),
	("-t", OptionF(str, 1, None, "Input topology PDB of GRO file")),
	("-stride", OptionF(int,1,1000, "Only read every Nth frame (default: 1000)")),
	("-om", OptionF(str, 1, "conformations.pdb", "Output PDB file with conformations"))
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

	optionsMotion = [
	# option for motion feature
	"Input/output options for project feature",
	("-p", OptionF(str, 1, None, "Input projection file")),
	("-aver", OptionF(str, 1, None, "PDB of GRO file for average structure")),
	("-ievec", OptionF(str, 1, None, "Eigenvector input file"))
	]

	optionsProjDist = [
	# options for projdist feature
	"Input/output options for projdist feature",
	("-p", OptionF(str, 1, None, "Input projection file")),
	("-pr", OptionF(str, 1, None, "Range of input projection files: \n\
	example: -pr proj_1.xvg-proj_10.xvg"))
	]

	optionsAutot = [
	# options for autot feature
	"Input/output options for autot feature",
	("-p", OptionF(str, 1, None, "Input projection file")),
	("-pr", OptionF(str, 1, None, "Range of input projection files: \n\
	example: -pr proj_1.xvg-proj_10.xvg")),
	("-ln", OptionF(int, 1, 1, "Number of lipids in the system (default: 1)")),
	("-dt", OptionF(float, 1, 0.01, "Timestep (ns) (default: 0.01 ns)")),
	("-o", OptionF(str, 1, "acor.xvg", "Name of output files"))
	]

	optionsKsst = [
	# options for ksst feature
	"Input/output options for ksst feature",
	("-p", OptionF(str, 1, None, "Input projection file")),
	("-pr", OptionF(str, 1, None, "Range of input projection files: \n\
	example: -pr proj_1.xvg-proj_10.xvg")),
	("-ln", OptionF(int, 1, 1, "Number of lipids in the system (default: 1)")),
	("-dt", OptionF(float, 1, 0.01, "Timestep (ns) (default: 0.01 ns)")),
	("-o", OptionF(str, 1, "kss.xvg", "Name of output files"))
	]

	optionsCombtrajs = [
	# options for combtrajs feature
	"Input/output options for combtrajs feature",
	("-fs", OptionF(str, -1, None, \
		"Input trajectories and corresponding average structures")),
	("-ou", OptionF(str,1,"united.xtc","Output combined trajectory")),
	("-oc", OptionF(str,1,"concatenated.xtc",\
		"Aligned trajectories for different simulations")),
	("-oa", OptionF(str,1,"average.pdb","Output average structures"))
	]

	options = [
	# list of all available features
	"List of procedures\n",
	"Performing PCA on lipid molecule conforamtions\n",
	("concat", Option(str,"con_traj",optionsConcat,\
		"Concatenate trajectories of individual lipids")),
	("conspace", Option(str,"lipic",optionsConspace,\
		"Vizualize possible conformations")),
	("covar", Option(str,"PCA",optionsCovar,\
		"Perform PCA on aligned concatenated lipid trajectory")),
	("project", Option(str,"project",optionsProject,\
		"Project concatenated lipid trajectory on the calculated PCs")),
	("motion", Option(str,"visualizing",optionsMotion,\
		"Create pdb file that represents the motion along selected PC")),
	("projdist", Option(str,"proj_dist_s",optionsProjDist,\
		"Plot projection distributions for selected PCs")),
	"\nCalculating characteristic timescales\n",
	("autot", Option(str,"autocorr",optionsAutot,\
		"Calculate autocorrelation decay times")),
	("ksst", Option(str,"kss_time_scale",optionsKsst,\
		"Calculate distribution convergence times")),
	"\nComparing several trajectories\n",
	("combtrajs", Option(str,"combtraj",optionsCombtrajs,\
		"Combine two trajectories"))
	]

	# add link for script; add link for publication
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
		if optionsF[ar].num == -1:
			listOfInputs = []
			while args:
				ar1 = args.pop(0)
				if ar1 in list(optionsF.keys()):
					optionsF[ar].setvalue(listOfInputs)
					args.insert(0,ar1)
					break
				else:
					listOfInputs.append(ar1)
			optionsF[ar].setvalue(listOfInputs)
		else:
			optionsF[ar].setvalue([args.pop(0) for i in range(optionsF[ar].num)]) # set value

	# Pass all the parameters to function
	params=list(v.value for v in optionsF.values())
	main(*params)




# 	elif args[0] == 'eigenvecdot':
# 		main = eigenvecdot.main
# 		if '-evec' in args:
# 			main(args[args.index('-evec') + 1], args[args.index('-evec') + 2])
# 		elif '-h' not in args and '-help' not in args:
# 			print('Missing parameters, try -h for flags\n')
# 		else:
# 			print('-evec <two files with eigenvectors (example: -evec eigenvec1.xvg eigenvec2.xvg)>')


# 	elif args[0] == 'pearson':
# 		main = pearson.main
# 		if '-cov1' in args and '-cov2' in args:
# 			main(args[args.index('-cov1') + 1], args[args.index('-cov2') + 1])
# 		elif '-h' not in args and '-help' not in args:
# 			print('Missing parameters, try -h for flags\n')
# 		else:
# 			print('-cov1, -cov2 - 2 files with covariance matrices')


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

if __name__ == '__main__':
	args = sys.argv[1:]
	main(args)
