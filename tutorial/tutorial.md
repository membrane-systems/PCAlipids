Below you could find the brief tutorial on using the software. The tutorial is organized as follows:
* Installation
* Analysing of single trajectory
* Comparing two different simulations

All the files needed to perform the tutorial could be downloaded [here](https://github.com/membrane-systems/PCAlipids/tree/master/tutorial/).

## Installation

### Prerequisites

PCALipids is a Python3 based software. To use it please install the following packages and programms:

* [Python 3.x](https://www.python.org/download/releases/3.0/)
* [Cython](http://cython.org/) - c-extensions for python (version >= 0.26.1)
* [Numpy](http://www.numpy.org/) - module for linear algebra (version >= 1.13.3)
* [Scipy](https://www.scipy.org/) - module for scientific calculations (version >= 0.19.1)
* [Matplotlib](https://matplotlib.org/) - module for creating plots (version >= )
* [MDTraj](http://mdtraj.org/1.9.0/) - module for MD trajectories routines (version >= 1.9.1)
* [Nose](http://nose.readthedocs.io/en/latest/) - module for python unittests (version >= 1.3.7)

The usefull information on installing Python and it's packages could be found here:
* https://wiki.python.org/moin/BeginnersGuide/Download
* https://packaging.python.org/tutorials/installing-packages/

We reccomend to install packages in the mentioned order to overcome possible issues.

### Welcome to PCALipids 

PCAlipids software is ready to use. You only need to download *scr* directory from the ripository. To run the software from any folder, you can add the path to the PCALipids on your machine to the PATH environmental variable:

    $ PATH=<path to PCALipids dir>:${PATH} 

where *<path to PCALipids dir>* has to be replaced with the path to PCALipids directory on your machine.

To run the software simply type in command prompt:

    $ pcalipids

The following output line should appear:

    $ Use -h or -help for more information

To get the information on usage type:

    $ pcalipids -h
    or
    $ pcalipids -help

You will get all the information about PCALipids functionality and usage. You could also find all the infromation [here](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

# Analysis of single trajectory

In this part of the tutorial we will work with the trajectories of DOPC lipid molecules. The prepared DOPC bilayer trajectory and structure could be found in the directory:

    $ tutorial/1_analysis_single/

Overall, the process is as follows:
* Trajectories of individual lipids are extracted from the input trajectory and are concatenated in a single trajectory
* The resulting trajectory is subjected to PCA. Covariance matrix, eigenvalues, eigenvectors and projections of the trajectory on eigenvectors are calculated.
* Projections of the trajectory on PCA eigenvectors are analyzed by calculating the projection distribution and characteristic time scales.

### Step 1: Creating concatenated trajectory

When working with lipids (or other flexible molecules) usually we have several molecules of interest in the system simulated (e.g. in the lipid bilayer there each lipid is a molecule of interest). To study available conformations it is worth to concatenate the trajectories of individual molecules into the concatenated trajectory. Thus, all of the possible conformations will be represented in this concatenated trajectory. The concatenated trajectory is produced by running:

    $ pcalipids concat -f trajectory.xtc -t structure.pdb -l DOPC

To get the information on *concat* procedure you can run

    $ pcalipids concat -h
    
or adress [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

After concatinating the trajectories of individual molecules the trajectory (*concatenated.xtc*) and the average structure (*average.pdb*) files are produced. To visualize possible conformations run

    $ pcalipids conspace -f concatenated.xtc -t average.pdb -stride 100
   
To get the information on *conspace* procedure you can run

    $ pcalipids conspace -h
    
or adress [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

As a result of *conspace* pruduces the *pdb* file with conformations of lipid molecule. The conformations could be visualized using PyMol or any graphical software you like. The example of lipid molecule conformations visualization is shown below.

![Example of comformational space with average structure](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/conform.png)

#### Step 2: Performing PCA

Now we are ready to move on to the next step. To perform PCA on the lipid molecule conformatoins run:

    $ pcalipids covar -f concatenated.xtc -t average.pdb
 
To get the information on *covar* procedure you can run

    $ pcalipids covar -h
    
or adress [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

This will calculate the covarience matrix, its eigenvectors and eigenvalues. You could visualize eigenvalues using any graphical package (matplotlib is great for python). You should get something similat to what you see below for eigenvalues and cumulitive eigenvalues:

![Eigenvalues of the covariance matrix](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/eigenvalues.png)
![Cumulative eigenvalues of the covariance matrix](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/eigenvalues_cumulative.png)

When the eigenvectors are calculated, we can project the trajectory on them:

    $ pcalipids project -f concatenated.xtc -t average.pdb -ia average.pdb -ievec eigenvec.xvg - first 1 -last 10 -op proj.xvg

To get the information on *project* procedure you can run

    $ pcalipids project -h
    
or adress [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

We have just calculated the trajectory projections on first 10 principal components. To visualize the motion along specific principal component (here we visualize the 1st PC) run:

    $ pcalipids motion -p proj.xvg -npc 1 -aver average.pdb -ievec eigenvec.xvg
    
To get the information on *motion* procedure you can run

    $ pcalipids motion -h
    
or adress [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

The *motion* produces 3 files:
* extreme_1_min.pdb - the structure projected on the 1st PC with the minimal projection value
* extreme_1_max.pdb - the structure projected on the 1st PC with the maximal projection value
* extreme1.pdb - 20 intermidiate structures representing the movement along the 1st PC

The movement along PC could be visualized using PyMol or any graphical software you like. The example of lipid molecule PC1 visualization is shown below.

![Example of single lipid motions](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/motion_PC_1.png)

The distribution of projections on principal components can be visualized using projdist:

    $ pcalipids projdist -p proj.xvg -first 1 -last 5

To get the information on *projdist* procedure you can run

    $ pcalipids projdist -h
    
or adress [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

*projdist* calculates the probability distribution functions for the projectoins of the trajectory on principal components. The resultiong PDFs are saved as *.png* files. See below the distributions for the first 2 PCs. 

![Probability distribution density of projections on the first principal component](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/PC1_dist.png)
![Probability distribution density of projections on the second principal component](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/PC2_dist.png)

**Covariance matrix** is a square matrix with the number of dimensions equal to the number of degrees of freedom in the system, and consequently divisible by 3. For output, the matrix is converted into one-dimensional array and written in blocks of 3 values per line.

**Eigenvalues** is a one-dimensional array of float numbers, each line contains a single eigenvalue.

**Eigenvectors** is a two-dimensional array of float numbers, each line contains a single eigenvector.

**Projections data**  is a one-dimensional array of float numbers for each of the components. The components are written one after another, each line contains the projection value for a particular projection for a particular frame.

#### Step 3: Characteristic timescales

To describe the equilibration process of the molecule of interest we could calculate the autocorrelation decay of the trajectory projections on the PC, or to study the process of PDFs convergence using Kolmogorov-Smirnov statistics. PCALipids can perform both types of the analysis. 

First, we need to split the projection file into the files, where only the projections on the specific PC will be:

    $ pcalipids splitproj -p proj.xvg
    
To get the information on *splitproj* procedure you can run

    $ pcalipids splitproj -h
    
or adress [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

This will produce 10 projection files *proj_i.xvg*. To calculate the autocorrelation decay times of the projections call *autot* procedure:

    $ pcalipids autot -pr proj_1.xvg-proj_10.xvg -ln 128 -dt 1

To get the information on *autot* procedure you can run

    $ pcalipids autot -h
    
or adress [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

This will produce 2 files and 2 png figures:
 * Computed autocorrelations for each principal component (*AUTO_VS_T.xvg* and *autocorrelation.png*)
 * Computed characteristic decay time for each principal component (*autocorr_relaxtime_vs_PC.xvg* and *autocorrelation_relax_time.png*)

You should get something similar to what you see below:

![Autocorrelation decay of the partial projections on differrent principal components](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/auto.png)
![Characteristic autocorrelation decay timescales](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/auto_TS.png)

To calculate the convergence of PDFs for the particular PC call *ksst* procedure:

    $ pcalipids ksst -pr proj_1.xvg-proj_10.xvg -ln 128 -dt 1

To get the information on *ksst* procedure you can run

    $ pcalipids ksst -h
    
or adress [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

This will produce 2 files and 2 png figures:
 * Computed KSS for each principal component (*KSS_vs_T.xvg* and *kss.png*)
 * Computed characteristic KSS convergence time for each principal component (*KSS_relaxtime_vs_PC.xvg* and *kss_relax_time.png*)

You should get something similar to what you see below:

![Examples of convergence of the distributions of the projections on principal components](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/KSS.png)
![Timescales](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/KSS_TS.png)

# Comparing two different simulations

In this part of the tutorial we will compare the simulations of DOPC lipid molecules conducted at different temperatures (310K & 338K). All the files for this work could be found in the directory
    
    $ tutorial/2_compare_simulations/

#### Step1: Preparing trajectories

First, we need to concatenate the trajectories for lipids in each simulation. You can use *concat* procedure for it, as done in the first part of the tutorial. It is important to save the resulting concatenated trajectories and average structures for each simulation in different files. Use optoins *-oc* and *-oa* to do it. To get more information on *concat* call

    $ pcalipids concat -h
    
or adress the [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt). You should get 2 different files with concatenated trajectories (*concatenated1.xtc* & *concatenated2.xtc*) and 2 different average structures (*average1.pdb* & *average2.pdb*).

To directly compare different simulations we need to compare them in the same basis. Thus, we first neen to concatenate different trajectories:

    $ pcalipids combtrajs -fs concatenated1.xtc average1.pdb concatenated2.xtc average2.pdb
    
To get more information of *combtrajs* call

    $ pcalipids combtrajs -h
    
or adress the [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

After combining the trajectories you should get the following files:
1. united.xtc - the united trajectory of two concatenated trajectories for different simulations
2. average.pdb - the average structure for the united trajectory
3. concatenated1_FALL.xtc - trajectory corresponding to the first simulation
4. average1_FALL.pdb - average structure for concatenated1_FALL.xtc
5. concatenated2_FALL.xtc - trajectory corresponding to the second simulation
6. average2_FALL.pdb - average structure for concatenated2_FALL.xtc

Now we are ready to perform PCA analysis on the united trajectory and compare the simulations.

#### Step2: Comparing the conformational spaces for different simulations

To compare the available conformations for different simulations we could apply different kinds of analysis:
* Direct comparison of covariance matrices
* Comparison of the sets of eigenvectors
* Comparison of the PDFs of the trajectory projections on the PCs

Two perform two first comparisons, we have to analyze the trajectories in their own basis. For the last comparison it is important to perform the analysis in the common basis.

##### Comparison of covariance matrices

To directly compare the covariance matrices we can use Pearson correlation coefficient. But first we need to calculate the covariance matrices for the concatenated trajectories of each simulation. You can do this by applying the *covar* procedure on *concatenated1_FALL.xtx* and *concatenated2_FALL.xtc* as it was done in the first part of the tutorial. Note, that you need to use respective average structure for the concatenated trajectory. Save the resulting covariance matrices, eigenvalues and eigenvectors to distinct files for each simulation (eg *cov1.dat*, *eigenval1.xvg* and *eigenvec1.xvg* for the first simulation; *cov2.dat*, *eigenval2.xvg* and *eigenvec2.xvg* for the second simulation). To get information on *procedure* call

    $ pcalipids covar -h
    
or adress the [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

Next, we can calculate the Pearson correlation coefficient of obtained covarience matrices:

    $ pcalipids pearson -cov1 cov1.dat -cov2 cov2.dat
    
To get more information on *pearson* procedure call

    $ pcalipids pearson -h
    
or adress the [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).
    
You should get the number in terminal close to **0.99916**. That means, that the available conformations of the lipid molecules at different temperatures are almost identical.

##### Comparison of sets of eigenvectors

To compare sets of eigenvectors for each simulation we can calculate the dot product matrix of two sets:

    $ pcalipids eigenvecdot -evec eigenvec1.xvg eigenvec2.xvg 

To get more information on *eigenvecdot* procedure call

    $ pcalipids eigenvecdot -h
    
or adress the [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

You should get the *eigenvecdot.dat* and *eigenvecdot.png* with the resulting dot product matrix. It should be similar to what you see below

![Scalar projections of evigenvectors from different trajectiories](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/1_vs_2/dot.png)

##### Comparison of simulations in common basis

To compare the PDFs of the trajectory projections on PCs we need to perform PCA on united trajectory. Use the *covar* procedure to do it. We get the covariance matrix (covar.dat), eigenvalues (eigenval.dat) and eigenvalues (eigenvec.dat). Then we need to project each concatenated trajectory on the common basis. For this, we can use *project* procedure applied to respective concatenated trajectories. You can find more on the *project* procedure in the first part of the tutorial. Note, that it is important to provide *project* with the eigenvectors obtained for the united trajectory. Finaly, we obtain a projection files for each concatenated trajectory (*proj1.xvg* and *proj2.xvg*). You can chose the PCs for which you want to perform the following analysis using options (*-first* and *-last*). We consider only first 10 modes for the following analysis.

Now we can plot the PDFs for the first PC for different simulations. First, we need to split each projection file into single projeciton files using *splitproj* procedure like we did in the first part of the tutorial. Then we can plot the PDF of the projections on the first PC

    $ pcalipids projdistm -file1 proj1_1.xvg -file2 proj2_1.xvg
    
To get more information on *projdistm* call

    $ pcalipids projdistm -h
    
or adress the [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

You should get the PNG figure with 2 PDFs for different simulations similar to what you see below

![PDFs for PC1 for different simulations analysed in common basis](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/1_vs_2/Distributions_in_general_basis.png)

#### Step3: Comparing timescales

To compare the timescales we will use procedures *autot* and *ksst* like we did in the first part of the tutorial. To plot the characteristic timescales for different simulations at the same plot use 

    $ pcalipids timescalespic -file1 proj1_autocorr_relaxtime_vs_PC_proj1.xvg -file2 proj2_autocorr_relaxtime_vs_PC_proj1.xvg -type auto -t t2
    
    $ pcalipids timescalespic -file1 proj1_KSS_relaxation_time_vs_PC.xvg -file2 proj2_KSS_relaxation_time_vs_PC.xvg -type kss -t t2

You shoul get the figures similar to what you see below

![Autocorrelation decay times for different simulations analysed in common basis](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/1_vs_2/Autocorrelation_relaxation_2_traj.png)

![KSS convergence times for different simulations analysed in common basis](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/1_vs_2/KSS_relaxation_2_traj.png)

To directly compare the timescales we can use *reltime* procedure

    $ pcalipids reltime -evec eigenval.xvg -time1 proj1_autocorr_relaxtime_vs_PC.xvg -time2 proj2_autocorr_relaxtime_vs_PC.xvg
    
    $ pcalipids reltime -evec eigenval.xvg -time1 proj1_KSS_relaxation_time_vs_PC.xvg -time2 proj2_KSS_relaxation_time_vs_PC.xvg
    
For autocorrelations you should get something close to **0.3374** and for KSS **0.59468**. Thus at higher temperatures the dynamics is speeded up by **~1.7**.
