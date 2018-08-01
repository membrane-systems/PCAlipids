1. 'concat' - create concatenated trajectory

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

2. 'covar' - principal component analysis

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

3. 'project' - calculating projections

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

4. 'projdist' - probability density

-p <projection file> (file format *.xvg)
 -first <first projection> -last <last projection> (int format). 
If not supplied, the first 3 projections will be analyzed.

5. 'ksst' - Kolmogorov-Smirnov convergence

-p <sequence of projection files> - this param must be the first
 -pr <range of files: "proj1.xvg-proj100.xvg">
-o <timescales file> (*.xvg)
-ln <number of lipids>
 -dt <timestep in (ns)>

6. 'autot' - Autocorrelation decay

-p <sequence of projection files> - this param must be the first
 -pr <range of files: "proj1.xvg-proj100.xvg">
-o <timescales file> (*.xvg)
-ln <number of lipids>

7. 'eigenvecdot' - scalar product of eigenvectors from different trajectories

8. 'conspace' - conformational space of lipids in trajectory

-f <trajectory file> (file format *.xtc)
-t <topology file> (any file with topology)
 -stride <<positive integer; step of reading frames>
 -om <output file with conformations>

9. 'pearson' - Pearson coefficient for covariance matrices from different trajectories

10. 'splitproj' - split files with all projections into files with projections for single PC

-p - projection file

11. 'motion' - demonstrates the motion along PC

-p - projection file
 -npc - number of principal component
 -aver - average structure
 -e - file with eigenvectors

