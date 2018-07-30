# PCALipids

Principal components analysis a standart tool to study the conformational changes in molecule kinetics. Usually it is applied for the systems with several metastable states (local minimas of free energy landscapes) and helps to develop the possible transitions between this states and corresponding transition rates. But it also can be usefull to study the conformations of flexible molecules (the free energy lanscape has only 1 sharp local minima). PCAlipids is a python based software that enables to perform deep quantitative analysis of conformations and dynamics of flexible molecules. All the information about the approach could be found in the following papers:

* [Principal Component Analysis of Lipid Molecule Conformational Changes in Molecular Dynamics Simulations, Buslaev et al., JCTC 2016](doi.org/10.1021/acs.jctc.5b01106)
* [Effects of Coarse Graining and Saturation of Hydrocarbon Chains on Structure and Dynamics of Simulated Lipid Molecules, Buslaev & Gushchin, Sci. Rep. 2017](doi.org/10.1038/s41598-017-11761-5)

Below you could find the brief tutorial on using the software. The tutorial is organized as follows:
* Installation
* Analysing of single trajectory
* Comparing two different simulations

## Installation

### Prerequisites

PCALipids is a Python3 based software. To use it please install the following packages and programms:

* [Python 3.x](https://www.python.org/download/releases/3.0/)
* [Cython]
* [Numpy](http://www.numpy.org/) - module for linear algebra (version >= 1.13.3)
* [Scipy](https://www.scipy.org/) - module for scientific calculations (version >= 0.19.1)
* [Matplotlib](https://matplotlib.org/) - module for creating plots
* [MDTraj](http://mdtraj.org/1.9.0/) - module for MD trajectories routines (version >= 1.9.1)
* [Nose]

The usefull information on installing Python and it's packages could be found here:
* https://wiki.python.org/moin/BeginnersGuide/Download
* https://packaging.python.org/tutorials/installing-packages/

We reccomend to install packages in the mentioned order to overcome possible issues.

### Welcome to PCALipids 

PCAlipids software is ready to use. You only need to download *scr* directory from the ripository. To run the software from any folder, you can add the path to the PCALipids on your machine to the PATH environmental variable:

   PATH=<path to PCALipids dir>:${PATH} 

where <path to PCALipids dir> has to be replaced with the path to PCALipids directory on your machine.

To run the software simply type in command prompt:

    $ pcalipids

The following output line should appear:

    $ Use -h or -help for more information

To get the information on usage type:

    $ pcalipids -h
    or
    $ pcalipids -help

You will get all the information about PCALipids functionality and usage. You could also find all the infromation here[https://github.com/membrane-systems/PCAlipids/blob/master/manual.txt]. **Need to add result of -h version to manual file**

**ANALYSIS**:

Please find below a step-by-step tutorial on using the software. Overall, the process is as follows:
* Trajectories of individual lipids are extracted from the input trajectory and are concatenated in a single trajectory
* The resulting trajectory is subjected to PCA. Covariance matrix, eigenvalues, eigenvectors and projections of the trajectory on eigenvectors are calculated.
* Projections of the trajectory on PCA eigenvectors are analyzed by calculating the projection distribution and characteristic relaxation times.

**NOTE**: Exemplary trajectory from the test_input directory can be used for tutorial.

#### Step 1: Creating concatenated trajectory

You need to place the PCAlipids script file in the folder that contains the trajectory (.xtc, .trr, etc.) and structure  (.pdb) files for your system. In our case, the names of the trajectory and structure files are “trajectory.xtc” and “structure.pdb”, respectively. The concatenated trajectory is produced by running:

    $ pcalipids concat -f trajectory.xtc -t structure.pdb

**Description**: Creates a concatenated trajectory.

**Input**: Trajectory file and structure file. Optional: reference structure for alignment, starting frame, frame step.

**Output**: Trajectory file and topology file for single lipid molecule.

**Parameters**:

**Required**:
* -f \<input trajectory file> 
* -t \<input topology file> 

**Optional**:
* -ref \<reference structure>
* -stride \<positive integer; step of reading frames> 
* -dt \<time in ps; number to determine from which frame to read the trajectory>
* -oc \<output trajectory file> - concatenated trajectory
* -oa \<output topology file> - average structure calculated from the concatenated trajectory
* -sf \<time in ps> - start frame for reading the trajectory
* -ef \<time in ps> - end frame for reading the trajectory


#### Step 2: Performing PCA

Now you are ready to move on to the next step. Run in the command prompt:
```bash 
$ pcalipids covar -f concatenated.xtc -t average.pdb
```
**Description**: Carry out the PCA of the concatenated trajectory.

**Input**: Concatenated trajectory file and structure file.

**Output**: Files with eigenvalues, eigenvectors and covariance matrix.

**Parameters**:

**Required**:
* -f \<input trajectory file> 
* -t \<input topology file>

**Optional**:
* -oeval \<output file with eigenvalues>
* -oevec \<output file with eigenvectors>
* -ocov \<output file with covariance matrix> 

**Files**:
* covar.dat - covariance matrix
* eigenval.xvg - eigenvalues
* eigenvec.xvg - eigenvectors

```bash 
$ pcalipids project -f concatenated.xtc -t average.pdb -ia ref_struct.pdb -ievec eigenvec.xvg - first 1 -last 10 - op proj.xvg
```
**Description**: Calculates projections.

**Input**: Concatenated trajectory file, structure file, average structure and eigenvectors. Optional: two positive integers to defining the range of principal components for the analysis.

**Output**: File with projections.

**Parameters**:

**Required**:
* -f \<input trajectory file> 
* -t \<input topology file>
* -ia \<input average structure>
* -ievec \<input eigenvectors>

**Optional**:
* -first \<number of the first principal component> 
* -last \<number of the last principal component>
* -op \<output file with projections>

**Files**:
* proj.xvg - projections of trajectory on principal components

We wrote all analyzing data in text files, so let us familiarize the structure of this files.

**Covariance matrix** is a square matrix with the number of dimensions equal to the number of degrees of freedom in the system, and consequently divisible by 3. For output, the matrix is converted into one-dimensional array and written in blocks of 3 values per line.

**Eigenvalues** is a one-dimensional array of float numbers, each line contains a single eigenvalue.

**Eigenvectors** is a two-dimensional array of float numbers, each line contains a single eigenvector.

**Projections data**  is a one-dimensional array of float numbers for each of the components. The components are written one after another, each line contains the projection value for a particular projection for a particular frame.

#### Step 3: Data processing for visualizing results

* The distribution of projections on principal components can be visualized using projdist:
```bash 
$ pcalipids projdist -p \<projection_file> -first \<number of the first projection> -last \<... last projection>
```
The output is a png plot with the distribution:
![Probability distribution density of projections on the first principal component](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/PC1_dist.png)
![Probability distribution density of projections on the second principal component](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/PC2_dist.png)

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

* Individual conformations can be visualized using "conspace" and "motion":
```bash 
$ pcalipids conspace -f \<concatenated trajectory file> -t \<average structure>
```
![Example of comformational space with average structure](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/1b.png)

* Single lipid motion along principal component can be vizualized using "motion":
```bash 
$ pcalipids motion -p \<projection file> -npc \<principal component> -aver \<average structure> -e \<file with eigenvectors>
```
![Example of single lipid motions](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/omegaasdasd.png)

* The equilibration of your lipid system can be easily examined by "ksst" and "autot" programs:
```bash 
$ pcalipids autot -pr proj_1.xvg-proj_100.xvg -ln 128 -dt 1
```
and
```bash    
$ pcalipids ksst -pr proj_1.xvg-proj_100.xvg -ln 128 -dt 1
```

**KSS-timesalces:**

![Examples of convergence of the distributions of the projections on principal components](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/KSS_tut.png "Examples of convergence of the distributions of the projections on principal components")
![Timescales](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/KSS_tut_relax_true.png "Characteristic distribution convergence timescales")

**Autrocorrelation-timescales:**

![Autocorrelation of the projections on differrent principal components](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/autot_tut.png "Autocorrelation of the projections on differrent principal components")
![Characteristic autocorrelation decay timescales](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/autot_tut_relax_1.png "Characteristic autocorrelation decay timescales")

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

