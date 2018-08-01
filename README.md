# PCALipids

Principal components analysis a standart tool to study the conformational changes in molecule kinetics. Usually it is applied for the systems with several metastable states (local minimas of free energy landscapes) and helps to develop the possible transitions between this states and corresponding transition rates. But it also can be usefull to study the conformations of flexible molecules (the free energy lanscape has only 1 sharp local minima). PCAlipids is a python based software that enables to perform deep quantitative analysis of conformations and dynamics of flexible molecules. All the information about the approach could be found in the following papers:

* [Principal Component Analysis of Lipid Molecule Conformational Changes in Molecular Dynamics Simulations, Buslaev et al., JCTC 2016](doi.org/10.1021/acs.jctc.5b01106)
* [Effects of Coarse Graining and Saturation of Hydrocarbon Chains on Structure and Dynamics of Simulated Lipid Molecules, Buslaev & Gushchin, Sci. Rep. 2017](doi.org/10.1038/s41598-017-11761-5)

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

where <path to PCALipids dir> has to be replaced with the path to PCALipids directory on your machine.

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
 * Computed autocorrelations for each principal component (*AUTO_VS_T.xvg* and **figure**)
 * Computed characteristic decay time for each principal component (*autocorr_relaxtime_vs_PC.xvg* and **figure**)

You should get something similar to what you see below:

![Autocorrelation decay of the partial projections on differrent principal components](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/auto.png)
![Characteristic autocorrelation decay timescales](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/auto_TS.png)

To calculate the convergence of PDFs for the particular PC call *ksst* procedure:

    $ pcalipids ksst -pr proj_1.xvg-proj_10.xvg -ln 128 -dt 1

To get the information on *ksst* procedure you can run

    $ pcalipids ksst -h
    
or adress [manual](https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt).

This will produce 2 files and 2 png figures:
 * Computed KSS for each principal component (*KSS_vs_T.xvg* and **figure**)
 * Computed characteristic KSS convergence time for each principal component (*KSS_relaxtime_vs_PC.xvg* and **figure**)

You should get something similar to what you see below:

![Examples of convergence of the distributions of the projections on principal components](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/KSS.png)
![Timescales](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/KSS_TS.png)

# Comparing two different simulations

* Two different trajectories can be compared using a single number - Pearson correlation coefficient of the respective covariance matrices - using programs "pearson" and "eigevecdot":
```bash 
$ pcalipids pearson -cov1 \<first file with covariance matrix> -cov2 \<second file with cov. matrix>
```
* The principal components obtained in different simulations can be compared using dot products of the respective eigenvectors:
```bash 
$ pcalipids eigenvecdot -evec \<first file with eigevector> \<second file>
```
The output plot is the dot product matrix (values in range (0; 1) -> (white; black)):
![Scalar projections of evigenvectors from different trajectiories](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/eigenveccomp.png)

## Contributing

Please read [CONTRIBUTING.txt](CONTRIBUTING.txt) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

Look for [VERSION.txt](VERSION.txt)

## Contacts

* **Khalid Mustafin** - khalid.mustafin@phystech.edu *developer of PCAlipids; code-related questions*
* **Pavel Buslaev** - pbuslaev@phystech.edu *applying PCAlipids to study lipids*
* **Ivan Gushchin** - ivan.gushchin@phystech.edu *general questions*

See also the list of [contributors](https://github.com/membrane-systems) who participated in this project.

## License

## Acknowledgments

