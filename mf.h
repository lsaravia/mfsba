/*  Copyright 2011 Leonardo A. Saravia
 
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
 
*/
#ifndef MF_H
#define MF_H
#include "smattpl.h"

int MultifractalSBA(simplmat <double> &pixval,simplmat <double> &q,char * outFile,
	int minBoxSize, int maxBoxSize, int deltaBoxSize=2,char normalize='S');

int MultifractalSBA(simplmat <double> &pixval,simplmat <double> &q,char * outFile,
	int minBoxSize, int maxBoxSize, int deltaBoxSize,char normalize, char * ident);

int MultifractalEBA(simplmat <double> &pixval,simplmat <double> &q,char * outFile,
	int minBoxSize, int maxBoxSize, int deltaBoxSize=2,char normalize='S');

int CoherenceLength(simplmat <double> &pixval,char * outFile,
	int minBoxSize, int maxBoxSize, int deltaBoxSize=2);

int PairCorrelation(simplmat <double> &pixval, char * outFile,
	double minBoxSize, double maxBoxSize, double deltaBoxSize, int maxPoints=0);

int MoranIRook(simplmat <double> &data, const char * outFile, const char * ident);
int PatchStats(simplmat <double> &data,int numSpecies, const char * outFile,const char * ident);

int MultispeciesReordering(simplmat <double> &data, simplmat <double> &newdata );

#endif  // MF_H

