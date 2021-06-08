Below you can find the tutorial on using PCAlipids. The tutorial is organized as follows:
* Installation
* Analysis of a single trajectory
* Comparison of two different trajectories

All the files needed to perform the tutorial can be downloaded [here](https://github.com/membrane-systems/PCAlipids/tree/master/tutorial/).

## Installation

### Prerequisites

PCAlipids is a Python3 based software. To use it please install the following packages and programs:

* [Python 3.x](https://www.python.org/download/releases/3.0/)
* [Cython](http://cython.org/) - c-extensions for Python (version >= 0.26.1)
* [Numpy](http://www.numpy.org/) - module for linear algebra (version >= 1.13.3)
* [Scipy](https://www.scipy.org/) - module for scientific calculations (version >= 0.19.1)
* [Matplotlib](https://matplotlib.org/) - module for creating plots (version >= 2.0.2)
* [MDTraj](http://mdtraj.org/1.9.0/) - module for MD trajectories routines (version >= 1.9.1)
* [Nose](http://nose.readthedocs.io/en/latest/) - module for Python unit tests (version >= 1.3.7)

Useful information on installing Python and its packages can be found here:
* https://wiki.python.org/moin/BeginnersGuide/Download
* https://packaging.python.org/tutorials/installing-packages/

We recommend to install packages in the mentioned order to overcome possible issues.

### Welcome to PCAlipids 

PCAlipids software is ready to use. You only need to download the *scr* directory from the repository. To run PCAlipids from any folder, you can add the path to the software on your machine to the PATH environmental variable:

    $ PATH=<path to PCALipids dir>:${PATH} 

where `<path to PCALipids dir>` has to be replaced with the path to PCAlipids directory on your machine.

To run the software simply type in the command prompt (on some machines you need to use pcalipids.py):

    $ pcalipids

As an output the description of the script functionality line should appear.

To get this information you can also type:

    $ pcalipids -h
    or
    $ pcalipids -help

You will get all the information about PCAlipids functionality and usage. You could also find all the information [here](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

# Analysis of a single trajectory

In this part of the tutorial we will work with a trajectory containing DOPC lipid molecules. The files (trajectory and structure) can be found in the directory:

    $ tutorial/1_analysis_single/

Overall, the process is as follows:
* Trajectories of individual lipids are extracted from the input trajectory and are concatenated in a single trajectory.
* The resulting trajectory is subjected to PCA. Covariance matrix, eigenvalues, eigenvectors, and projections of the trajectory on eigenvectors are calculated.
* Projections of the trajectory on PCA eigenvectors are analyzed by calculating the projection distribution and characteristic time scales.

### Step 1: Creating concatenated trajectory

When working with lipids (or other flexible molecules) usually we have several molecules of interest in the simulated system (e.g. in the lipid bilayer, where each lipid is a molecule of interest). To study the available conformations, it is worth to concatenate the trajectories of individual molecules into a concatenated trajectory. Thus, all the possible conformations will be represented in this concatenated trajectory. The concatenated trajectory is produced by running:

    $ pcalipids concat -f trajectory.xtc -t structure.pdb -l DOPC

To get the information on *concat* procedure you can run

    $ pcalipids concat -h
    
or address the [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

After concatenating the trajectories of individual molecules, the trajectory (*concatenated.xtc*) and the average structure (*average.pdb*) files are produced. To visualize possible conformations, run

    $ pcalipids conspace -f concatenated.xtc -t average.pdb -stride 100
   
To get the information on *conspace* procedure you can run

    $ pcalipids conspace -h
    
or address [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

As a result, *conspace* produces the *pdb* file with conformations of lipid molecule. The conformations can be visualized using PyMol or any graphical software you like. An example of lipid molecule conformations visualization is shown below.

![Example of comformational space with average structure](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/conform.png)

#### Step 2: Performing PCA

Now we are ready to move on to the next step. To perform PCA on lipid molecule conformations run:

    $ pcalipids covar -f concatenated.xtc -t average.pdb
 
To get the information on *covar* procedure you can run

    $ pcalipids covar -h
    
or address [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

This will calculate the covariance matrix, its eigenvectors and eigenvalues. You could visualize eigenvalues using any graphical package (matplotlib is great for Python). You should get something similar to what you see below for eigenvalues and cumulative eigenvalues:

![Eigenvalues of the covariance matrix](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/eigenvalues.png)
![Cumulative eigenvalues of the covariance matrix](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/eigenvalues_cumulative.png)

After the eigenvectors are calculated, we can project the trajectory on them:

    $ pcalipids project -f concatenated.xtc -t average.pdb -ia average.pdb -ievec eigenvec.xvg -first 1 -last 10 -op proj.xvg

To get the information on *project* procedure you can run

    $ pcalipids project -h
    
or address [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

We have just calculated the trajectory projections on the first 10 principal components. To visualize the motion along specific principal component (here we visualize the 1st PC) run:

    $ pcalipids motion -p proj_1.xvg -aver average.pdb -ievec eigenvec.xvg
    
To get the information on *motion* procedure you can run

    $ pcalipids motion -h
    
or address [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

The *motion* produces 1 file:
* extreme1.pdb - 20 intermediate structures representing the movement along the 1st PC

The movement along a PC can be visualized using PyMol or any graphical software you like. An example of such a visualization is shown below.

![Example of single lipid motions](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/motion_PC_1.png)

The distribution of projections on principal components can be visualized using projdist:

    $ pcalipids projdist -p proj_1.xvg 

or

    $ pcalipids projdist -pr proj_1.xvg-proj_3.xvg

To get the information on *projdist* procedure you can run

    $ pcalipids projdist -h
    
or address [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

*projdist* calculates the probability distribution functions (PDFs) for the projections of the trajectory on principal components. The resulting PDFs are saved as *.png* files. See below the distributions for the first 2 PCs. 

![Probability distribution density of projections on the first principal component](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/PC1_dist.png)
![Probability distribution density of projections on the second principal component](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/PC2_dist.png)

**Covariance matrix** is a square matrix with the number of dimensions equal to the number of degrees of freedom in the system, and consequently divisible by 3. For output, the matrix is converted into one-dimensional array and written in blocks of 3 values per line.

**Eigenvalues** is a one-dimensional array of float numbers, where each line contains a single eigenvalue.

**Eigenvectors** is a two-dimensional array of float numbers, where each line contains a single eigenvector.

**Projections data**  is a one-dimensional array of float numbers for each of the components. The components are written one after another, and each line contains the projection value for a particular projection for a particular frame.

#### Step 3: Characteristic timescales

To describe the equilibration process of the molecule of interest we either calculate the decay in autocorrelation of the PC projection values, or look at the convergence of PC projections' PDFs using Kolmogorov-Smirnov statistics (KSS). PCAlipids can perform both types of the analysis. 

To calculate the autocorrelation decay times of the projections run the *autot* procedure:

    $ pcalipids autot -pr proj_1.xvg-proj_10.xvg -ln 128 -dt 1

To get the information on *autot* procedure you can run

    $ pcalipids autot -h
    
or address [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

This will produce 2 files and 2 png figures:
 * Computed autocorrelations for each principal component (*AUTO_VS_T.xvg* and *autocorrelation.png*)
 * Computed characteristic decay time for each principal component (*autocorr_relaxtime_vs_PC.xvg* and *autocorrelation_relax_time.png*)

You should get something similar to what you see below:

![Autocorrelation decay of the partial projections on differrent principal components](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/acor_values_vs_t.png)
![Characteristic autocorrelation decay timescales](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/acor_relax_times_vs_pc.png)

To calculate the convergence of PDFs for a particular PC run *ksst* procedure:

    $ pcalipids ksst -pr proj_1.xvg-proj_10.xvg -ln 128 -dt 1

To get the information on *ksst* procedure you can run

    $ pcalipids ksst -h
    
or address [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

This will produce 2 files and 2 png figures:
 * Computed KSS for each principal component as a function of time (*KSS_vs_T.xvg* and *kss.png*)
 * Computed characteristic KSS convergence time for each principal component (*KSS_relaxtime_vs_PC.xvg* and *kss_relax_time.png*)

You should get something similar to what you see below:

![Examples of convergence of the distributions of the projections on principal components](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/kss_values_vs_t.png)
![Timescales](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/kss_relax_times_vs_pc.png)

# Comparing two different simulations

In this part of the tutorial we will compare two simulations of DOPC lipid molecules conducted at two different temperatures (310 K & 338 K). All the files for this work can be found in the directory
    
    $ tutorial/2_compare_simulations/

#### Step1: Preparing the trajectories

First, we need to concatenate the trajectories for lipids in each simulation. You can use *concat* procedure to do that, similarly to the first part of the tutorial. It is important to save the resulting concatenated trajectories and average structures for each simulation in different files. Use options *-oc* and *-oa* to do it. To get more information on *concat* run

    $ pcalipids concat -h
    
or address the [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt). You should get 2 different files with concatenated trajectories (*concatenated1.xtc* & *concatenated2.xtc*) and 2 different average structures (*average1.pdb* & *average2.pdb*).

To directly compare different simulations we need to compare them in the same PC basis. Thus, we need first to concatenate different trajectories:

    $ pcalipids combtrajs -fs concatenated1.xtc average1.pdb concatenated2.xtc average2.pdb
    
To get more information of *combtrajs* run

    $ pcalipids combtrajs -h
    
or address the [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

After combining the trajectories you should get the following files:
1. united.xtc - the united trajectory of two concatenated trajectories for different simulations
2. average.pdb - the average structure for the united trajectory

Now we are ready to perform PCA analysis on the concatenated trajectory and to compare the simulations.

#### Step2: Comparing the conformational spaces for different simulations

To compare the available conformations for different simulations we can apply different kinds of analysis:
* Direct comparison of covariance matrices
* Comparison of the sets of eigenvectors
* Comparison of the PDFs of the trajectory projections on the PCs

To perform the first two comparisons, we have to analyze the trajectories in their own basis. For the last comparison it is important to perform the analysis in the common basis.

##### Comparison of covariance matrices

To directly compare the covariance matrices we can use the Pearson correlation coefficient. But first we need to calculate the covariance matrices for the concatenated trajectories of each simulation. This can be done by applying the *covar* procedure on *concatenated1.xtx* and *concatenated2.xtc* as it was done in the first part of the tutorial. Note, however, the following:
* Your two trajectories must have the same order of atoms. (This might not be the case, if you are comparing trajectories from two different force fields.)
* You should preferably use the respective average structure for each of your two concatenated trajectories. (Otherwise the mean atom positions in your trajectories will not be zero.)
* These two average structures must be aligned to one another, because otherwise the entries in the two covariance matrices might not correspond to the same coordinates. (So you need to do two extra alignment steps: first to align the two average structures, and second to realign the corresponding trajectories to the newly aligned average structures using the switch *-r* in *concat*.)
Save the resulting covariance matrices, eigenvalues and eigenvectors to different files for each simulation (eg *cov1.dat*, *eigenval1.xvg* and *eigenvec1.xvg* for the first simulation; *cov2.dat*, *eigenval2.xvg* and *eigenvec2.xvg* for the second simulation). To get information on *covar* run

    $ pcalipids covar -h
    
or address the [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

Next, we can calculate the Pearson correlation coefficient of the obtained covarience matrices:

    $ pcalipids pearson -cov1 cov1.dat -cov2 cov2.dat
    
To get more information on *pearson* procedure run

    $ pcalipids pearson -h
    
or address the [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).
    
The resulting value as reported in the terminal should be close to **0.99916**. This means that the conformations available to the lipid molecules at different temperatures are almost identical.

##### Comparison of sets of eigenvectors

To compare sets of eigenvectors for each simulation we can calculate the dot product matrix of the two sets:

    $ pcalipids evecdot -evec eigenvec1.xvg eigenvec2.xvg 

To get more information on *eigenvecdot* procedure run

    $ pcalipids evecdot -h
    
or address the [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

You should get the *eigenvecdot.dat* and *eigenvecdot.png* with the resulting dot product matrix. It should be similar to what you see below:

![Scalar projections of evigenvectors from different trajectiories](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/1_vs_2/dot.png)

##### Comparison of simulations in common PC basis

To compare the PDFs of the trajectory projections on common PCs, we need to perform PCA on the united concatenated trajectory that contains the lipids of both simulations (the *united.xtc* you obtained using *combtrajs* in Step1 above). Use the *covar* procedure on this trajectory to perform the PCA. We get the covariance matrix (*covar.dat*), eigenvalues (*eigenval.dat*), and eigenvalues (*eigenvec.dat*).

Then we need to project the concatenated trajectories of each simulation (the *concatenated1.xtc* & *concatenated2.xtc* you obtained in Step1 above) on the common basis. For this, we can use the *project* procedure applied to the respective concatenated trajectories. You can find more on the *project* procedure in the first part of the tutorial. Note that it is important to provide *project* with the eigenvectors obtained for the united trajectory. Finally, we obtain projection files for each concatenated trajectory (*proj1.xvg* and *proj2.xvg*). You can choose the PCs for which you want to perform the following analysis by using the options (*-first* and *-last*); we consider only the first 10 modes for the following analysis.

Now we can plot the PDFs for the first PC for different simulations. To plot the PDF of the projections on the first PC run

    $ pcalipids projdistm -files proj1_1.xvg proj2_1.xvg
    
To get more information on *projdistm* run

    $ pcalipids projdistm -h
    
or address the [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

You should get a PNG figure with 2 PDFs for different simulations similar to what you see below:

![PDFs for PC1 for different simulations analysed in common basis](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/distributions.png)

#### Step3: Comparing timescales

To compare the timescales we will use procedures *autot* and *ksst* like we did in the first part of the tutorial. You may specify name of output files for different trajectories (e.g. *acor1.xvg*, *acor2.xvg*, *kss1.xvg*, *kss2.xvg*). To plot the characteristic timescales for different simulations at the same plot use 

    $ pcalipids tsCmpFig -file1 acor1_relaxtime_vs_pc.xvg -file2 acor2_relaxtime_vs_pc.xvg -type auto -t t2
    
    $ pcalipids tsCmpFig -file1 kss1_relaxation_time_vs_pc.xvg -file2 kss2_relaxation_time_vs_pc.xvg -type kss -t t2

You should get the figures similar to what you see below:

![Autocorrelation decay times for different simulations analysed in common basis](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/1_vs_2/Autocorrelation_relaxation_2_traj.png)

![KSS convergence times for different simulations analysed in common basis](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/1_vs_2/KSS_relaxation_2_traj.png)

To directly compare the timescales we can use *reltime* procedure. It is recommended to use relaxation time for all PCs, however the script will work with the specified range of projections as well (in case of this tutorial - first 10 PCs)

    $ pcalipids reltime -eval eigenval.xvg -time1 acor1_relaxtime_vs_pc.xvg -time2 acor2_relaxtime_vs_pc.xvg
    
    $ pcalipids reltime -eval eigenval.xvg -time1 kss1_relaxtime_vs_pc.xvg -time2 kss2_relaxtime_vs_pc.xvg
    
For autocorrelations you should get something close to **0.358** and for KSS **0.558**. Thus, at higher temperatures the dynamics is accelerated by **~1.7**.
