/*  Copyright 2012 Leonardo A. Saravia
 
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

#include <math.h>
#include <ctype.h>
#include <string>
#include "RWFile.h"
#include "mf.h"
using namespace std;


int MultispeciesReordering(simplmat <double> &data, simplmat <double> &newdata )
{
	int i,dimY,dimX,numSpecies=0,a,maxi;
	dimX = data.getRows();
	dimY = data.getCols();
	if( newdata.getRows()!=dimX || newdata.getCols()!=dimY)
		newdata.resize(dimX,dimY,0.0);
	else
		newdata.fill(0.0);
		
	double totCells=dimX*dimY,maxDen;
	
	for(i=0; i<dimY; i++)
		for(int j=0;j<dimX;j++)
			{
			a = data(j,i);
			if( a>numSpecies )
				numSpecies=a;
			}
	if( numSpecies<1 || numSpecies>totCells)
		return(0);
		
	double * den = new double[numSpecies];

	for(i=0; i<numSpecies; i++)
		den[i]=0;
	
	for(i=0; i<dimY; i++)
		for(int j=0;j<dimX;j++)
			{
			a = data(j,i);
			if( a>0 )
				den[ a-1 ]++;
			}
			
	int newSpecie=0;
	while(newSpecie<numSpecies)
	{
		maxDen=maxi=0;
		for( i=0; i<numSpecies; i++)
			if(den[i]>maxDen)
				{
					maxi=i+1;
					maxDen=den[i];
				}
		if (maxDen==0) break;
		
		cout << maxi << "-" << maxDen << "\t";
		newSpecie++;
		for(i=0; i<dimY; i++)
			for(int j=0;j<dimX;j++)
				{
				a = data(j,i);
				if( a==maxi )
					newdata(j,i)=newSpecie;
				}
		den[maxi-1]=0;
	}
	return(1);
}

