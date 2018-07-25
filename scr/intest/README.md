# PCALipids

LIMOS is a program for analyzing lipid trajectories. Global properties of lipid bilayer such as area per lipid, bilayer thickness, density and etc do not provide information about structural properties of lipid molecules. LIMOS is based on the method called principal component analysis, which can be easily applied to trajectories with one lipid molecule. As a result we can get a quantitative description of the conformational space of a single lipid molecule. 

## Info

The current tutorial is the step-by-step instruction for using LIMOS program package. The exercises that are described here aim to help you understand the basics of data analysis applied to MD trajectories of lipids. The main idea is to provide the approach to the study of lipids dynamics. It goes on to describe the pipeline of these analysis, starting from converting trajectories of lipid bilayers into trajectory with single lipid molecule, finishing with some information on the interpretation of the data. So, the time that you will spend on this tutorial depends on the computational facilities and your own interest and sense of purpose. Good luck and have fun!

### Prerequisites

Before doing this tutorial you should install Python interpreter on your working machine. Also if you wish to do the examples in this tutorial, you must have some python modules installed on your computer. Below are the modules that you need to install for the correct work of LIMOS, you can also click on to find tutorials on installing these packages and more useful information:

* [Interpreter of Python 3.x](https://www.python.org/download/releases/3.0/)
* [Numpy](http://www.numpy.org/) - helpful module for linear algebra
* [Scipy](https://www.scipy.org/) - module for scientific calculations
* [Mdtraj](http://mdtraj.org/1.9.0/) - reading, writing and analyzing MD trajectories 

If you have already installed the above modules on your computer we suggest you go directly to the LIMOS tutorial. 

### Installing

It is advisable to install the modules in the correct order, so you will avoid the subsequent difficulties. 

First of all you should install Python interpreter version 3.x. If you are using Ubuntu 14.10 or newer, then you  can install Python (recommended 3.5+) by the following commands in command prompt:

    $ sudo apt update
    $ sudo apt install python3

To see which version of Python you have installed run next command:

    $ python3 --version

So you have successfully installed Python 3 interpreter, to launch it run

    $ python3

or to execute python script run

    $ python3 script.py

**NOTE**: How install third-party modules packages and modules for python?
The most crucial third-party Python package is pip. Once installed, you can download and install any modules for Python. Python 3.4 and later versions include pip by default, so to check if pip installed, open command prompt and run:
    $ command -v pip 

It is important to note that pip package belongs to Python interpreter version 2.x and pip3 package belongs to Python 3.x. 


### Some PCALipids basics:

To run LIMOS on your computer open command prompt and run:

    $ python3 limos.py

Obviously, nothing happens, except that an output line appears:

    $ Use -h or -help for more information

Let’s follow the advice and enter the following command in prompt:

    $ python3 limos.py -h
    or
    $ python3 limos.py -help

This command displays information about the functions that are implemented in the limos.
Let's try to analyze the trajectory with LIMOS!

**Step by step**:

Analyzing the lipid bilayer trajectory with LIMOS consists basically of 2 principal steps:
* Converting the trajectory of lipid bilayer into trajectory with single lipid molecule without water
* Providing principal component analysis for concatenated trajectory

#### Step 1

You need to place the LIMOS script file in the folder that contains the trajectory (*.xtc, *.trr, etc.) and file of the topology (*.pdb) of your system.
In our case, name of trajectory file - “trajectory.xtc” and topology file - “topology.pdb”. To start, as it was said, it is necessary to create a concatenated trajectory. Run next command:

    $ python3 limos.py concat -f trajectory.xtc -t topology.pdb

We specified the trajectory file and the topology file for our script, but did not specify the names of the output files, so they will be called “concatenated.xtc” and “average.pdb” for trajectory and topology, respectively. Also, in your working folder, files of the concatenated trajectory and topology should appear.

#### Step 2

Now you are ready to move on to the next step. Run in the command prompt:

    $ python3 limos.py covar -f concatenated.xtc -t average.pdb

After execution you will find 4 new files in your working directory.

Files:
* covar.dat - covariance matrix
* eigenval.xvg - eigenvalues
* eigenvec.xvg - eigenvectors
* projection.xvg - projections of trajectory on principal components

We wrote all analyzing data in text files, so let us familiarize the structure of this files.

**Covariance matrix** has equal dimensions and they always are divisible by 3 (because of degrees of freedom are divisible by 3 too), so we wrote it in file using the following algorithm:
    1) Transformation the covariance matrix into one-dimensional array
    2) Divide the array into blocks of length 3
    3) Write each block inline

**Eigenvalues** data is a one-dimensional array of float numbers, so we write them line by line in file.

**Eigenvectors** data is a two-dimensional array of float numbers, so we write to the line only one component, one after the other.

**Projections data**  is a one-dimensional array of float numbers, the number of projections equals to number of frames in concatenated trajectory. We write them one after one in file.

### Data processing for visualizing results


## Scientific base

[]

## Contributing

Please read [CONTRIBUTING.txt](https://github.com/KhalidMustafin/ScAns/blob/master/limos/CONTRIBUTING.txt) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

Look for [VERSION.txt](https://github.com/KhalidMustafin/ScAns/blob/master/limos/VERSION.txt)

## Authors

* **Ivan Gushchin** - *Laboratory head, initial idea*
* **Pavel Buslaev** - *Developer, initial idea*
* **Khalid Mustafin** - *Developer*

See also the list of [contributors](https://github.com/membrane-systems) who participated in this project.

## License

This project is licensed ...

## Acknowledgments

* ...