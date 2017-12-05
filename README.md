# PCAlipids

This is a software for analyzing lipid trajectories using PCA. The analysis results in comprehensive description of conformations and dynamics of lipid molecules. The methodology is based on following papers:
* [Principal Component Analysis of Lipid Molecule Conformational Changes in Molecular Dynamics Simulations, Buslaev et al., JCTC 2016](doi.org/10.1021/acs.jctc.5b01106)
* [Effects of Coarse Graining and Saturation of Hydrocarbon Chains on Structure and Dynamics of Simulated Lipid Molecules, Buslaev & Gushchin, Sci. Rep. 2017](doi.org/10.1038/s41598-017-11761-5)

The software uses the approaches and terminology similar to GROMACS covar and anaeig utilities.

## Info

Please find below a step-by-step tutorial on using the software. Overall, the process is as follows:
* Trajectories of individual lipids are extracted from the input trajectory and are concatenated in a single trajectory
* The resulting trajectory is subjected to PCA. Covariance matrix, eigenvalues, eigenvectors and projections of the trajectory on eigenvectors are calculated.
* Projections of the trajectory on PCA eigenvectors are analyzed by calculating the projection distribution and characteristic relaxation times.

**NOTE**: Exemplary trajectory from the test_input directory can be used for tutorial.

### Prerequisites

The software requires a Python interpreter and some other modules that you need to install:

* [Interpreter of Python 3.x](https://www.python.org/download/releases/3.0/)
* [Numpy](http://www.numpy.org/) - helpful module for linear algebra
* [Scipy](https://www.scipy.org/) - module for scientific calculations
* [MDTraj](http://mdtraj.org/1.9.0/) - reading, writing and analyzing MD trajectories 

### Installing Python

It is advisable to install the modules in the correct order, so you will avoid subsequent difficulties. 

First of all you need to install the Python interpreter version 3.x. If you are using Ubuntu 14.10 or newer, then you can install Python (recommended 3.5+) using the following commands in command prompt:

    $ sudo apt update
    $ sudo apt install python3

To see which version of Python you have installed run next command:

    $ python3 --version

To launch the Python 3 interpreter run:

    $ python3

To execute python script (script.py) run:

    $ python3 script.py

**NOTE**: Third-party python modules and module packages can be downloaded and installed using *pip*. Python 3.4 and later versions include pip by default, so to check if pip installed, open command prompt and run:
    
    $ pip -V 

It is important that *pip* package belongs to Python interpreter version 3.x. If this is not the case, please install python3-pip. Then you can go to the downloaded PCAlipids directory and run:

    $ sudo pip install -r requirements.txt

### Installing PCAlipids

PCAlipids does not require installation. You can download the files and use them as is. To run the software from any folder, add the path to the software directory to the PATH global variable in ~/.bashrc as follows:

   PATH=/path/to/your/pcalipids/directory:${PATH} 

where /path/to/your/pcalipids/directory has to be replaced with your path to your pcalipids directory

### PCAlipids basics:

To run the software on your computer open command prompt and run:

    $ pcalipids

The following output line will appear:

    $ Use -h or -help for more information

Let’s follow the advice and enter the following command in prompt:

    $ pcalipids -h
    or
    $ pcalipids -help

This command displays information about the functions that are implemented in the PCAlipids.
Let's try to analyze the trajectory using it!

**ANALYSIS**:

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

#### Step 2: Performing PCA

Now you are ready to move on to the next step. Run in the command prompt:
```bash 
    $ pcalipids covar -f concatenated.xtc -t average.pdb
```
**Description**: Carry out the PCA of the concatenated trajectory.

**Input**: Concatenated trajectory file and structure file. Optional: two positive integers to defining the range of principal components for the analysis.

**Output**: Files with eigenvalues, eigenvectors, covariance matrix and projections.

**Parameters**:

**Required**:
* -f \<input trajectory file> 
* -t \<input topology file>

**Optional**:
* -first \<number of the first principal component> 
* -last \<number of the last principal component>
* -oeval \<output file with eigenvalues>
* -oevec \<output file with eigenvectors>
* -ocov \<output file with covariance matrix> 
* -op \<output file with projections>

**Files**:
* covar.dat - covariance matrix
* eigenval.xvg - eigenvalues
* eigenvec.xvg - eigenvectors
* projection.xvg - projections of trajectory on principal components

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
The output plot is the dot product matrix (values in range (-1; +1) -> (dark blue; yellow)):
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
$ pcalipids autot -pr proj_1.xvg-proj_100.xvg -ln 128 -dt 0.01
```
and
```bash    
$ pcalipids ksst -pr proj_1.xvg-proj_100.xvg -ln 128 -dt 0.01
```
![Examples of convergence of the distributions of the projections on principal components](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/KSS_tut.png "Examples of convergence of the distributions of the projections on principal components")
![Timescales](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/KSS_tut_relax.png "Characteristic distribution convergence timescales")

![Autocorrelation of the projections on differrent principal components](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/autot_tut.png "Autocorrelation of the projections on differrent principal components")
![Characteristic autocorrelation decay timescales](https://github.com/membrane-systems/PCAlipids/blob/master/scr/output/autot_tut_relax.png "Characteristic autocorrelation decay timescales")

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

