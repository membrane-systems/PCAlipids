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
import overlap

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
	("-om", OptionF(str, 1, "conformations.pdb", "Output PDB file with conformations")),
	("-al", OptionF(int, 1, 0, "Is align needed (0: no; 1: yes; default: 0)"))
	]

	optionsCovar = [
	# options for covar feature
	"Input/output options for covar feature",
	("-f", OptionF(str, 1, None, "Input XTC or TRR concatenated trajectory")),
	("-t", OptionF(str, 1, None, "Input topology PDB of GRO file")),
	("-oeval", OptionF(str, 1, "eigenval.xvg", "Eigenvalue output file")),
	("-oevec", OptionF(str, 1, "eigenvec.xvg", "Eigenvector output file")),
	("-ocov", OptionF(str, 1, "covar.dat", "Covariance matrix output file"))
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
	"Input/output options for motion feature",
	("-p", OptionF(str, 1, None, "Input projection file")),
	("-aver", OptionF(str, 1, None, "PDB of GRO file for average structure")),
	("-ievec", OptionF(str, 1, None, "Eigenvector input file")),
	("-eo", OptionF(str, 1, "extreme.pdb", \
		"PDB file name with a motion for selected PC"))
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

	optionsPearson = [
	# options for pearson feature
	"Input/output options for pearson feature",
	("-cov1", OptionF(str,1,None,"Covariance matrix for the first trajectory")),
	("-cov2", OptionF(str,1,None,"Covariance matrix for the second trajectory")),
	("-o", OptionF(str,1,"pearson.dat","Output file")),
	"Pearson correlation coeficient is shown in the terminal and saved to the output file"
	]

	optionsOverlap = [
	# options for pearson feature
	"Input/output options for overlap feature",
	("-cov1", OptionF(str,1,None,"Covariance matrix for the first trajectory")),
	("-cov2", OptionF(str,1,None,"Covariance matrix for the second trajectory")),
	("-o", OptionF(str,1,"overlap.dat","Output file")),
	"Overlap of the conformational spaces as introduce in https://doi.org/10.1103/PhysRevE.65.031910 \n\
	is shown in the terminal and saved to the output file"
	]

	optionsEigenvecdot = [
	# options for eigenvecdot feature
	"Input/output options for evecdot feature",
	("-evec", OptionF(str,2,None,"Input eigenvectors for 2 simulations")),
	("-o", OptionF(str,1,"eigenvecproduct.xvg","Output file name"))
	]

	optionsProjdistM = [
	# options for projdistm feature
	"Input/output options for projdistm feature",
	("-files", OptionF(str,-1,None,"Input projection files")),
	("-o", OptionF(str,1,"distributions.png","Output file name"))
	]

	optionsTimescalespic = [
	# options for timescalespic feature
	"Input/output options for timescalespic feature",
	("-file1", OptionF(str,-1,None,"Input timescales for 1st simulation")),
	("-file2", OptionF(str,-1,None,"Input timescales for 2nd simulation")),
	("-type", OptionF(str,1,"auto","auto for autocorrelation timescales;\n\
	kss for distribution convergence timescales")),
	("-t",  OptionF(str,-1,"t2",'t1" or "t2" for selected \n\
	timescale measure decay in e or e^2 times')),
	("-o", OptionF(str,1,"timescales_comp.png","Output file name"))
	]

	optionsReltime = [
	# options for reltime feature
	"Input/output options for reltime feature",
	("-eval",OptionF(str,1,None,"Eigenvalue file for base simulation")),
	("-time1",OptionF(str,1,None,"characteristic timescales for 1st simulation")),
	("-time2",OptionF(str,1,None,"characteristic timescales for 2nd simulation")),
	("-o",OptionF(str,1,"tsCmp.dat","Output file"))
	]

	optionsEvals = [
	# options for evals feature
	"Input/output options for evals feature",
	("-ieval", OptionF(str,1,None,"Eigenvalue file for simulation")),
	("-cum", OptionF(int,1,0,"Plot cumulative (or just) eigenvalues 1(0). default: 0")),
	("-o",OptionF(str,1,"eval.png","Output file"))
	]

	options = [
	# list of all available features
	"List of procedures\n",
	"Performing PCA on lipid molecule conforamtions\n",
	("concat", Option(str,"con_traj",optionsConcat,\
		"Concatenate trajectories of individual lipids")),
	("conspace", Option(str,"lipic",optionsConspace,\
		"Vizualize possible conformations for concatenated trajectory")),
	("covar", Option(str,"PCA",optionsCovar,\
		"Perform PCA on aligned concatenated lipid trajectory")),
	("evals", Option(str,"eigenvalues",optionsEvals,\
		"Plot eigenvalues for calculated PCs")),
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
		"Combine two trajectories")),
	("pearson", Option(str,"pearson",optionsPearson,\
		"Compare covariance matrices from two simulations")),
	("overlap", Option(str,"overlap",optionsOverlap,\
		"Compare covariance matrices from two simulations")),
	("evecdot", Option(str,"eigenvecdot",optionsEigenvecdot,\
		"Compare eigenvectors from two simulations")),
	("projdistm", Option(str,"projdistm",optionsProjdistM,\
		"Plot projection distributions for simulations of interest")),
	("tsCmpFig", Option(str,"timescalespic",optionsTimescalespic,\
		"Plot timescales for two trajectories")),
	("reltime", Option(str,"reltime",optionsReltime,\
		"Compare characteristic timescales for two trajectories"))
	]

	# Short description of the software
	desc = """
   PCAlipids is a software for analysis of lipid 
   molecule conformations and dynamics.
   You could find more information on the software usage at:
      https://github.com/membrane-systems/PCAlipids
   In case of usage for your research please cite:
   1. [Principal Component Analysis of Lipid Molecule Conformational 
   Changes in Molecular Dynamics Simulations, Buslaev et al., JCTC 2016]
   2. [Effects of Coarse Graining and Saturation of Hydrocarbon 
   Chains on Structure and Dynamics of Simulated Lipid Molecules, 
   Buslaev & Gushchin, Sci. Rep. 2017]
	"""

	# If the user asks for help: pcalipids.py -h
	if (len(args)>0 and (args[0] == '-h' or args[0] == '--help')) or len(args)==0:
		print("\n",__file__[__file__.rfind('/')+1:]) # print executable file name
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
	funcname=args[0]
	optionsF = options[args[0]].optList

	# We don't need 0 argument of args any more
	args = args[1:]

	# if the user asks for help: pcalipids.py concat -h
	if (len(args)>0 and (args[0] == '-h' or args[0] == '--help')) or len(args)==0:
		print("\n",__file__[__file__.rfind('/')+1:],\
			funcname+": "+options[funcname].description) # print executable file name
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

if __name__ == '__main__':
	args = sys.argv[1:]
	main(args)
