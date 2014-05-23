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

struct outRegress {

    double bt,ba,bf,rt2,ra2,rf2,sdbt,sdba,sdbf;
    outRegress():bt(0),ba(0),bf(0),rt2(0),ra2(0),rf2(0),sdbt(0),sdba(0),sdbf(0){};
};

int loglogRegress(simplmat <double> &q,int &numBoxSizes, simplmat <double> &box,
    simplmat <double> &tauQ, simplmat <double> &alphaQ, simplmat <double> &fQ, simplmat <outRegress> &oR);

int standardBoxCount(simplmat <double> &pixval,simplmat <double> &q, int &minBoxSize, 
    int &maxBoxSize, int &numBoxSizes, char &normalize,
    simplmat <double> &box, simplmat <double> &tauQ, simplmat <double> &alphaQ, simplmat <double> &fQ,
    double (*winMov)(simplmat <double> &pixval,const int &rowIni,const int &rowEnd,const int &colIni,const int &colEnd)
);

int standardBoxCountSAD(simplmat <double> &pixval,simplmat <double> &q, int &minBoxSize, 
    int &maxBoxSize, int &numBoxSizes, char &normalize,
    simplmat <double> &box, simplmat <double> &tauQ, simplmat <double> &alphaQ, simplmat <double> &fQ
);

double winMovSum(simplmat <double> &pixval,const int &rowIni,const int &rowEnd,const int &colIni,const int &colEnd);
double winMovNumSp(simplmat <double> &px,const int &rowIni,const int &rowEnd,const int &colIni,const int &colEnd);
double winMovSAD(simplmat <double> &px,const int &rowIni,const int &rowEnd,const int &colIni,const int &colEnd,const double &qT);

#endif  // MF_H

