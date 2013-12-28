//
//
#include "Randomizations.h"
#include <iomanip>
//#include "fortify.h"


void Randomizations::Randomize(simplmat<double>& data,simplmat<double>& rdata)
{
	int dimX = data.getRows();
	int dimY = data.getCols();
	rdata.resize(dimX,dimY);
	
	int max = dimX*dimY;
	int x,y,p=0;
	randomizePosXY * pos = new randomizePosXY[max];
	for(x=0; x<dimX; x++)
		for(y=0; y<dimY; y++)
		{
			pos[p].x=x;
			pos[p].y=y;
			p++;
		}

	
    random_shuffle(pos, pos + max, rnd);

//	JumbleXY(pos, max);
	p=0;
	for(x=0; x<dimX; x++)
		for(y=0; y<dimY; y++)
		{
			int px = pos[p].x;
			int py = pos[p].y;
			rdata(x,y) = data(px,py);
			p++;
		}
	
	delete []pos;
}


void Randomizations::Randomize(simplmat<double>& data)
{
	double * pdata = data.pointer();
	int max = data.getRows()*data.getCols();
    random_shuffle(pdata, pdata + max, rnd);
}
	


