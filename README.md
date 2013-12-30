
Multifractal estimation using a standard box-counting algorithm
===============================================================

Software used to determine multifractal spectra: Generalized dimensions Dq and spectrum 
of singularities f(alfa), the software was used in the paper:

﻿Saravia LA, Giorgi A, Momo F (2012) Multifractal growth in periphyton communities. Oikos 121: 1810–1820.

*mfSBA* estimate both spectra and outputs data to evaluate the fits. A more detailed
description is in the file mfSBA_README.

*mfSBArnz* estimate the Dq spectrun and randomize the original image N times to 
calculate a confidence interval to test the hipothesis that the original distribution
was random. 

*multiSpeciesSBA* estimates the multifractal spectra of a 2D species distribution assuming that each position is one individual and that each value represents a different specie.

*sed2grad* transforms a sed file using a discrete gradient transformation. This program was used in the paper: 

Saravia LA, Giorgi A, Momo FR (2012) Multifractal Spatial Patterns and Diversity in an Ecological Succession. PLoS One 7: e34096. 

The programs accept as input tiff, and sed files (ASCII format).

b4-991008bio.sed/tif  are examples from the paper of 2D biomass distributions to analyze with the program. 

q21.sed is a file with q values used to estimate spectra.

Please leave an Issue on github if you have some trouble.

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
 
