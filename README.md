
Multifractal estimation using a standard box-counting algorithm
===============================================================

Software used to determine multifractal spectra: Generalized dimensions Dq and spectrum 
of singularities f(alfa), the software is described in the paper:

1. Saravia LA (2014) mfSBA: Multifractal analysis of spatial patterns in ecological communities [v2; ref status: indexed, http://f1000r.es/347]. F1000Research 3: 14. doi:10.12688/f1000research.3-14.v2.

*mfSBA* estimate both spectra and outputs data to evaluate the fits. A more detailed
description is in the file mfSBA_README.

*mfSBArnz* estimate the Dq spectrun and randomize the original image N times to 
calculate a confidence interval to test the hipothesis that the original distribution
was random. 

*multiSpeciesSBA* estimates the multifractal spectra of a 2D species distribution assuming that each position is one individual and that each value represents a different specie. The species can be analyzed as a spatial rank surface (SRS) or as a generalization of species area distribution which I call multifractal SAD (species abundance distribution). 

*sed2grad* transforms a sed file using a discrete gradient transformation. 

### R Functions and scripts 

There are two scripts to test and demonstrate the software:

**testMFA.r**: for standard multifractal analysis


**testMFA_SAD.r**: for SRS and SAD multifractal analysis

Several ".sed" files ".tif" are used as examples

The **q21.sed** is a file with q values used to estimate spectra.



### These programs were used in the following papers: 

1. Saravia LA, Giorgi A, Momo FR (2012) Multifractal Spatial Patterns and Diversity in an Ecological Succession. PLoS One 7: e34096. 
﻿
2. Saravia LA, Giorgi A, Momo F (2012) Multifractal growth in periphyton communities. Oikos 121: 1810–1820.




> Please leave an Issue on github if you have some trouble.



License
=======

	Copyright 2011 Leonardo A. Saravia
 
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
