
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
   + A: Normalize measure and save
        

### Output:

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
have to divide it by (q-1) see equation 2 in [2]. The 3rd,4rd are the estimated alfa and f(alfa). Following are the coeficients of determination of the linear regresions for each of the previous, and finally the sd's.


# multispeciesSBA

The program assumes that the input file is a distribution of species, calculates the rank abundance distribution (RAD) and replace each species by its rank, this is called the species rank surface (SRS see [3]) then it applies multifractal spectra estimation using the canonical method [1]. 
The program also calculates the multifractal distribution of species abundances distribution (SAD), as a generalization of species area curve, this is described in detail in [4]. NOTE: for method E (Based on SAD) log(Zq) in t.[inputFile] and Tau in s.[inputFile] are multiplied by -1 so to calculate Dq=Tau/(q-1) like in the other methods.


### Usage:

	mf inputFile qFile minBox maxBox numBoxSizes option

The parameters and output files are the same as previous, the only difference is options

**option**
   + N: Not normalize measure use SRS
   + S: Normalize measure
   + A: Save the SRS distribution
   + E: Use the spatial SAD to compute dimensions
        This option implies N

 


#### To do:

* 1 dimension and quasi 1 dimension estimation is not implemented but it should be very easy to add it



### References 

1. Chhabra, Meneveau, Jensen & Sreenivasan, Phys.Rev.A 40,5284 (1989)

2. Saravia LA, Giorgi A, Momo F (2012) Multifractal growth in periphyton communities. Oikos 121: 1810–1820. doi:10.1111/j.1600-0706.2011.20423.x.

3. Saravia LA (2014) mfSBA: Multifractal analysis of spatial patterns in ecological communities [v2; ref status: indexed, http://f1000r.es/347]. F1000Research 3: 14. doi:10.12688/f1000research.3-14.v2.

4. Borda-de-Água L, Hubbell SP, McAllister M (2002) Species-Area Curves, Diversity Indices, and Species Abundance Distributions: A Multifractal Analysis. Am Nat 159: 138–155.