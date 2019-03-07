   PCAlipids is a software for analysis of lipid 
   molecule conformations and dynamics.
   You could find more information on the software usage at:
      https://github.com/membrane-systems/PCAlipids
   In case of usage for your research please cite:
   1. [Principal Component Analysis of Lipid Molecule Conformational 
   Changes in Molecular Dynamics Simulations, Buslaev et al., JCTC 2016](https://doi.org/10.1021/acs.jctc.5b01106)
   2. [Effects of Coarse Graining and Saturation of Hydrocarbon 
   Chains on Structure and Dynamics of Simulated Lipid Molecules, 
   Buslaev & Gushchin, Sci. Rep. 2017](https://doi.org/10.1038/s41598-017-11761-5)
	
List of procedures

Performing PCA on lipid molecule conforamtions

    concat: Concatenate trajectories of individual lipids
  conspace: Vizualize possible conformations
     covar: Perform PCA on aligned concatenated lipid trajectory
     evals: Plot eigenvalues for calculated PCs
   project: Project concatenated lipid trajectory on the calculated PCs
    motion: Create pdb file that represents the motion along selected PC
  projdist: Plot projection distributions for selected PCs

Calculating characteristic timescales

     autot: Calculate autocorrelation decay times
      ksst: Calculate distribution convergence times

Comparing several trajectories

 combtrajs: Combine two trajectories
   pearson: Compare covariance matrices from two simulations
   evecdot: Compare eigenvectors from two simulations
 projdistm: Plot projection distributions for simulations of interest
  tsCmpFig: Plot timescales for two trajectories
   reltime: Compare characteristic timescales for two trajectories

Detailed functions descriptions:

    concat: Concatenate trajectories of individual lipids

Input/output options for concat feature
        -f: Input trajectory file (.xtc, .trr, ...)
        -t: Input topology file (.pdb, .gro, ...)
   -stride: Only read every Nth frame
       -sf: First frame (ps) to read from trajectory
       -ef: Last frame (ps) to read from trajectory
       -oc: Output concatenated trajectory file
       -oa: Output average structure
        -r: Input reference file (.pdb, .gro). If not supplied, the structure of the first lipid and the first frame is used for alignment
        -l: Lipid type

  conspace: Vizualize possible conformations

Input/output options for conspace feature
        -f: Input XTC or TRR file
        -t: Input topology PDB of GRO file
   -stride: Only read every Nth frame (default: 1000)
       -om: Output PDB file with conformations

     covar: Perform PCA on aligned concatenated lipid trajectory

Input/output options for covar feature
        -f: Input XTC or TRR concatenated trajectory
        -t: Input topology PDB of GRO file
    -oeval: Eigenvalue output file
    -oevec: Eigenvector output file
     -ocov: Covariance matrix output file
-invertPC1: Invert the distribution for the first PC or not
		[default] 0 - do not invert
		1 - invert

   	 evals: Plot eigenvalues for calculated PCs

Input/output options for evals feature
    -ieval: Eigenvalue file for simulation
      -cum: Plot cumulative (or just) eigenvalues 1(0). default: 0
        -o: Output file

   project: Project concatenated lipid trajectory on the calculated PCs

Input/output options for project feature
        -f: Input XTC or TRR concatenated trajectory
        -t: Input topology PDB of GRO file
       -ia: PDB of GRO file for average structure
    -ievec: Eigenvector input file
    -first: First PC for projection (default: 1)
     -last: First PC for projection (default: 10)
       -op: Output projection files name
		Do not use '-' or '.' symbols in the projection file names

    motion: Create pdb file that represents the motion along selected PC

Input/output options for motion feature
        -p: Input projection file
     -aver: PDB of GRO file for average structure
    -ievec: Eigenvector input file

  projdist: Plot projection distributions for selected PCs
  
Input/output options for projdist feature
        -p: Input projection file
       -pr: Range of input projection files: 
	example: -pr proj_1.xvg-proj_10.xvg

     autot: Calculate autocorrelation decay times

Input/output options for autot feature
        -p: Input projection file
       -pr: Range of input projection files: 
	example: -pr proj_1.xvg-proj_10.xvg
       -ln: Number of lipids in the system (default: 1)
       -dt: Timestep (ns) (default: 0.01 ns)
        -o: Name of output files

      ksst: Calculate distribution convergence times

Input/output options for ksst feature
        -p: Input projection file
       -pr: Range of input projection files: 
	example: -pr proj_1.xvg-proj_10.xvg
       -ln: Number of lipids in the system (default: 1)
       -dt: Timestep (ns) (default: 0.01 ns)
        -o: Name of output files

 combtrajs: Combine two trajectories

Input/output options for combtrajs feature
       -fs: Input trajectories and corresponding average structures
       -ou: Output combined trajectory
       -oc: Aligned trajectories for different simulations
       -oa: Output average structures

   pearson: Compare covariance matrices from two simulations

Input/output options for pearson feature
     -cov1: Covariance matrix for the first trajectory
     -cov2: Covariance matrix for the second trajectory
        -o: Output file
Pearson correlation coeficient is shown in the terminal and saved to the output file

   evecdot: Compare eigenvectors from two simulations

Input/output options for pearson feature
     -evec: Input eigenvectors for 2 simulations
        -o: Output file name

 projdistm: Plot projection distributions for simulations of interest
 
Input/output options for evecdot feature
    -files: Input projection files
        -o: Output file name

  tsCmpFig: Plot timescales for two trajectories

Input/output options for tsCmpFig feature
    -file1: Input timescales for 1st simulation
    -file2: Input timescales for 2nd simulation
     -type: auto for autocorrelation timescales;
	kss for distribution convergence timescales
        -t: t1" or "t2" for selected 
	timescale measure decay in e or e^2 times
        -o: Output file name

   reltime: Compare characteristic timescales for two trajectories

Input/output options for reltime feature
     -eval: Eigenvalue file for base simulation
    -time1: characteristic timescales for 1st simulation
    -time2: characteristic timescales for 2nd simulation
        -o: Output file