#include "smattpl.h"
#include "RWFile.h"
#include <math.h>

int main(int argc, char * argv[])
{
	if( argc<2 )
	{	
		cerr << "Converts the file with a gradient transformation " << endl;
		cerr << "using the followint formula:" << endl;
		cerr << "sqrt( ((bi1,j−bi−1,j))^2 + (b(i,j1)−b(i,j−1))^2 )" << endl;
		cerr << "Usage: sed2grad input output" << endl;
		exit(1);
	}
	simplmat<double> data;
	RWFile file;
	file.ReadSeed(argv[1], data);
	int dimX = data.getRows();
	int dimY = data.getCols();
	int i,j;	
	// Convierte a Gradiente
	//
	//
	int grad = 1; // atoi(argv[3]);
	
	if( grad)
	{
		simplmat<double> ndat(dimX,dimY);
		double xx,yy;
		for(i=0; i<dimX; i++)
			for( j=0; j<dimY; j++)
			{
				if(i==dimX-1)
					xx = data(i, j) - data(i-1,j);
				else if(i==0)
				  xx = data(i+1, j) - data(i,j);
				else
					xx = data(i+1, j) - data(i-1,j);
					
				if(j==dimY-1)
					yy = data(i, j) - data(i,j-1);
				else if(j==0)
					yy = data(i, j+1) - data(i,j);
				else
					yy = data(i, j+1) - data(i,j-1);
				
				ndat(i,j)=sqrt( xx*xx + yy*yy );
			}
		for(i=0; i<dimX; i++)
			for( j=0; j<dimY; j++)
				data(i,j)=ndat(i,j);

   	}
	file.WriteSeed(argv[2],data);
}
