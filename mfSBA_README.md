
# mfSBA

Multifractal spectra estimation using the canonical method [1]

### Usage:

    mfSBA inputFile qFile minBox maxBox numBoxSizes option

**inputFile**: this file can have only two formats 1) one layer tiff 2) sed file format.
Sed File is a format I invented to use with my stochastic cellular automata models to
represent a square grid of values. 
it has a header of two lines, the first line describes the dimensions X Y of the data
the second line describes the type of data, for this program the type must be BI, that
means the grid can have real numbers with double precision. See the example files.

**qFile**: is a sed file with the q to calculate the multifractal spectra

**minBox,maxBox,numBoxSizes**: Minimun box size, maximun box size and number of sizes of boxes.
The program uses boxsize in powers of two, if maxBoxSize is greater than the half of the image size is set to that value.
numBoxSize: if is less than the number of powers of 2 that number will be used cutting from the bigger ones, so put a bigger number 

**option**
   + N: Not normalize measure
   + S: Normalize measure -> SUM all the pixels and divide each pixel by that value
   + D: Add 1 and normalize -> Add 1 to all the image then normalize
   + A: Normalize measure and save
        

Output:
-------
t.[inputFile]
The first line are labels, box size, log box size and q requested in [qFile]  
The first two columns are box size and Log box size used in estimations and the following columns are
log(Zq) (equation 1 in [2]) used to calculate Dq.

a.[inputFile]
Same as t file but for alpha.

f.[inputFile]
Same as t file but for f(alpha).

s.[inputFile]
The first line are labels.
The first column are q's requested, 2nd is Tau the slope of log(zq) vs log(box size) to calculate Dq you
have to divide it by (q-1) see equation 2 in [2]. The 3rd,4rd are the estimated alfa and f(alfa). Following 
are the coeficients of determination of the linear regresions for each of the previous, and finally the sd's.



To do:
Add MLE following 
﻿Newman, M. E. J. 2005. Power laws, Pareto distributions and Zipf’s law. Contemporary Physics 46:323-351.



[1]Chhabra, Meneveau, Jensen & Sreenivasan, Phys.Rev.A 40,5284 (1989)
[2]Saravia LA, Giorgi A, Momo FR (2012) Multifractal growth in periphyton communities. Oikos 0: 0.
