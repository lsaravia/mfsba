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

#include <math.h>
#include <ctype.h>
#include <string>
#include "RWFile.h"
#include "mf.h"
using namespace std;

// 	Q: Valores de Q posibles 
//	pixval: Measure
//  fileout: Archivo de Output
//
int MultifractalSBA(simplmat <double> &pixval,simplmat <double> &q, char * outFile,
	int minBoxSize, int maxBoxSize, int numBoxSizes, char normalize)
{
//	simplmat <double> piQ;
//	simplmat <double> piHat;
	simplmat <double> box;
	simplmat <double> tauQ;
	simplmat <double> alphaQ;
	simplmat <double> fQ;
	simplmat <int> boxIni;
	simplmat <int> boxFin;
	
	double sumAlphaQ, sumFQ,cnt,qT,piQT,piHatT,tauQT;

	int  boxSize, iq, i, iRow, iCol, ix, iy, qNum, xDim, yDim;
	int  yResto=0,xResto=0,actBoxSize=0;

	qNum = q.getRows();
	xDim = pixval.getRows();
	yDim = pixval.getCols();

	cnt=0;
   	double tot=0;

	for(iy=0; iy<yDim; iy++)
		for(ix=0; ix<xDim; ix++)
			cnt+=pixval(ix,iy);

    if(toupper(normalize)=='S')
    {
		for(iy=0; iy<yDim; iy++)
			for(ix=0; ix<xDim; ix++)
			{
				pixval(ix,iy)/=cnt;
				tot+=pixval(ix,iy);
			}
	}		
	else if(toupper(normalize)=='D')
	{
		cnt+=xDim*yDim;
		for(iy=0; iy<yDim; iy++)
			for(ix=0; ix<xDim; ix++)
			{
				pixval(ix,iy)=(pixval(ix,iy)+1)/cnt;					
				tot+=pixval(ix,iy);
			}
	}

	
	if( maxBoxSize > yDim/2 || maxBoxSize > xDim/2  )
		maxBoxSize = (xDim<yDim ? xDim : yDim)/2;

//
//	Calculo de boxes con intevalo constante en escala logaritmica
//
/*	double deltaBoxSize = (log10(static_cast<double>(maxBoxSize)) - log10(static_cast<double>(minBoxSize)))/ numBoxSizes;
	box.resize(numBoxSizes);

	box(0)= minBoxSize;
	for(i=1; i<numBoxSizes; i++)
	{
		box(i)=pow10(log10(box(i-1))+deltaBoxSize);

		if(box(i) > maxBoxSize )
		{
			numBoxSizes = i;
			break;
		}
		
	}
	for(i=0; i<numBoxSizes; i++)
		box(i) = ceil(box(i));
		
	for(i=numBoxSizes-1; i>1; i--)
	{
		if( box(i) == box(i-1) )
		{
			for(int ii=i; ii<numBoxSizes-1; ii++)
				box(ii) = box(ii+1);
			numBoxSizes--;
		}
	}
*/	
	// Calculo de boxes como potencias de 2
	//
	box.resize(numBoxSizes);
   	int finDimX=0,finDimY=0;

	for(finDimX=0,i=0;finDimX<minBoxSize;finDimX=pow(2.0,i))
    	i++;
    int minPot=i;
    
	box(0)= pow(2.0,minPot);
	for(i=1; i<numBoxSizes; i++)
	{
		box(i)=pow(2.0,minPot+i);
		if( box(i)>maxBoxSize )
		{
			numBoxSizes = i;
			break;
		}
	}

	for(finDimX=0;finDimX<xDim;finDimX=pow(2.0,i))
    	i++;
	finDimX=pow(2.0,i-1);
	for(finDimY=0,i=1;finDimY<yDim;finDimY=pow(2.0,i))
    	i++;
	finDimY=pow(2.0,i-1);

	int numRep;

	alphaQ.resize(numBoxSizes,qNum,0.0);
	tauQ.resize(numBoxSizes,qNum,0.0);
	fQ.resize(numBoxSizes,qNum,0.0);
	boxIni.resize(4,2,0);
	boxFin.resize(4,2,0);
	
	for(boxSize=0; boxSize<numBoxSizes; boxSize++)
	{
		for(iq=0;iq<qNum; iq++)
		{
			qT = q(iq);

			actBoxSize = box(boxSize);
			yResto = yDim - finDimY;
			xResto = xDim - finDimX;

			boxIni.fill(0);
			boxFin.fill(finDimX);
			boxFin(0,1)=boxFin(1,1)=boxFin(2,1)=boxFin(3,1)=finDimY;
			numRep=1;
			
			if( xResto>0 && yResto>0  )
			{
				numRep=4;
				boxIni(1,0)= xResto;
				boxIni(3,0)= xResto;
				boxIni(2,1)= yResto;
				boxIni(3,1)= yResto;
				boxFin(1,0)= finDimX+xResto;
				boxFin(3,0)= finDimX+xResto;
				boxFin(2,1)= finDimY+yResto;
				boxFin(3,1)= finDimY+yResto;

			}
			if( xResto>0  && yResto==0)
			{
				numRep=2;
				boxIni(1,0)= xResto;
				boxFin(1,0)= finDimX+xResto;
			}
			if( yResto>0  && xResto==0)
			{
				numRep=2;
				boxIni(1,1)= yResto;
				boxFin(1,1)= finDimY+yResto;
			}

			for(int rep=0;rep<numRep;rep++)
			{
				piQT=0;
				tauQT=0;
				
				for(iRow=boxIni(rep,1); iRow <= boxFin(rep,1)-actBoxSize; iRow+=actBoxSize )
					for(iCol=boxIni(rep,0); iCol <= boxFin(rep,0)-actBoxSize; iCol+=actBoxSize )
					{
						cnt = 0.0;

//						window examination
						for(iy=iRow; iy<iRow+actBoxSize; iy++)
							for(ix=iCol; ix<iCol+actBoxSize; ix++)
								cnt+=pixval(ix,iy);
						if( cnt>0.0 )
							piQT+=pow(cnt,qT);
						
					}
				
				sumAlphaQ=0;
				sumFQ=0;
				piHatT=0;

				if( piQT > 0.0 )
					tauQT=log10(piQT);
				else
				{
					goto AFTERSUM;
				}
			
				// To do AlphaQ and FQ

				//	window movement
				for(iRow=boxIni(rep,1); iRow <= boxFin(rep,1)-actBoxSize; iRow +=actBoxSize )
					for(iCol=boxIni(rep,0); iCol <=boxFin(rep,0)-actBoxSize; iCol +=actBoxSize )
					{

						cnt=0.0;

						//	window examination
						for(iy=iRow; iy<iRow+actBoxSize; iy++)
							for(ix=iCol; ix<iCol+actBoxSize; ix++)
								cnt+=pixval(ix,iy);

						if( cnt>0.0 )
						{
							piHatT=pow(cnt,qT)/piQT;
							sumAlphaQ+=piHatT*log10(cnt);
							sumFQ+=piHatT*log10(piHatT);
						}
						
//						if( (piQT != 0.0) && (cnt!=0.0) )
//							piHatT=pow(cnt,qT)/piQT;
//						else
//							piHatT=0;
//
//						if( cnt > 0.0 )
//							sumAlphaQ+=piHatT*log10(cnt);
//						
//						if(piHatT > 0.0)
//							sumFQ+=piHatT*log10(piHatT);
					}
				
				AFTERSUM:
				alphaQ(boxSize,iq)+=sumAlphaQ;
				fQ(boxSize,iq)+=sumFQ;
				tauQ(boxSize,iq)+=tauQT;
			}
			
			tauQ(boxSize,iq)/=static_cast<double>(numRep);
			alphaQ(boxSize,iq)/=static_cast<double>(numRep);
			fQ(boxSize,iq)/=static_cast<double>(numRep);
			
		}
	}


	string tFileName("t.");
	tFileName+=outFile;
	ofstream tFile(tFileName.c_str());
   
	string aFileName("a.");
	aFileName+= outFile;
	ofstream aFile(aFileName.c_str());
   
	string fFileName("f.");
	fFileName+= outFile;
	ofstream fFile(fFileName.c_str());

	tFile << "Box Size" << "\t" << "Log Box";
	aFile << "Box Size" << "\t" << "Log Box";
	fFile << "Box Size" << "\t" << "Log Box";
	for(i=0;i<qNum;i++)
	{
		tFile << "\t" << q(i);
		aFile << "\t" << q(i);
		fFile << "\t" << q(i);
	}
	tFile << endl;
	aFile << endl;
	fFile << endl;
	for(boxSize=0; boxSize<numBoxSizes; boxSize++)
	{
		tFile << box(boxSize) << "\t" << log10(box(boxSize));
		for(i=0;i<qNum;i++)
			tFile << "\t" << tauQ(boxSize,i);
		tFile << endl;
		
		aFile << box(boxSize) << "\t" << log10(box(boxSize));
		for(i=0;i<qNum;i++)
			aFile << "\t" << alphaQ(boxSize,i);
		aFile << endl;
		
		fFile << box(boxSize) << "\t" << log10(box(boxSize));
		for(i=0;i<qNum;i++)
			fFile << "\t" << fQ(boxSize,i);
		fFile << endl;
	}
	
	// Regression of the log - log variables
	//
	//
	string sFileName("s.");
	sFileName+=outFile;
	ofstream sFile(sFileName.c_str());
	sFile << "q\tTau\talfa\tf(alfa)\tR-Tau\tR-alfa\tR-f\tSD-Tau\tSD-alfa\tSD-f" << endl;
	for(i=0; i<qNum; i++)
	{
		double sumx=0,
			sumyt = 0,
			sumya = 0,
			sumyf = 0,
			sumxyt = 0,
			sumxya = 0,
			sumxyf = 0,
			sumxsq = 0,
			sumysqt = 0,
			sumysqa = 0,
			sumysqf = 0,
			sdbt=0,
			sdba=0,
			sdbf=0;
		
		for(boxSize=0; boxSize<numBoxSizes; boxSize++)
		{
			double x = log10(box(boxSize));
			double yt = tauQ(boxSize,i);
			double ya = alphaQ(boxSize,i);
			double yf = fQ(boxSize,i);
			
			sumx += x;
			sumyt += yt;
			sumya += ya;
			sumyf += yf;
			sumxyt += x * yt;
			sumxya += x * ya;
			sumxyf += x * yf;
			sumxsq += x * x;
			sumysqt += yt * yt;
			sumysqa += ya * ya;
			sumysqf += yf * yf;
		}
		double xbar  = sumx/numBoxSizes ;
		double ybart = sumyt/numBoxSizes ;
		double ybara = sumya/numBoxSizes ;
		double ybarf = sumyf/numBoxSizes ;
	
		double bt = (sumxyt - numBoxSizes *xbar*ybart)/(sumxsq - numBoxSizes *xbar*xbar);
		double ba = (sumxya - numBoxSizes *xbar*ybara)/(sumxsq - numBoxSizes *xbar*xbar);
		double bf = (sumxyf - numBoxSizes *xbar*ybarf)/(sumxsq - numBoxSizes *xbar*xbar);
	//	a = ybar - (b*xbar);
		double rt = (sumxyt - numBoxSizes *xbar*ybart)/sqrt((sumxsq - numBoxSizes  * xbar * xbar) * (sumysqt - numBoxSizes *ybart*ybart));
		double ra = (sumxya - numBoxSizes *xbar*ybara)/sqrt((sumxsq - numBoxSizes  * xbar * xbar) * (sumysqa - numBoxSizes *ybara*ybara));
		double rf = (sumxyf - numBoxSizes *xbar*ybarf)/sqrt((sumxsq - numBoxSizes  * xbar * xbar) * (sumysqf - numBoxSizes *ybarf*ybarf));
		if( numBoxSizes  > 4 )
		{
			sdbt = sqrt((1-rt*rt)*(sumysqt-numBoxSizes*ybart*ybart))/((numBoxSizes -4)*(sumxsq-numBoxSizes *xbar*xbar));
			sdba = sqrt((1-ra*ra)*(sumysqa-numBoxSizes*ybara*ybara))/((numBoxSizes -4)*(sumxsq-numBoxSizes *xbar*xbar));
			sdbf = sqrt((1-rf*rf)*(sumysqf-numBoxSizes*ybarf*ybarf))/((numBoxSizes -4)*(sumxsq-numBoxSizes *xbar*xbar));
		}
		sFile << q(i) << "\t" << bt << "\t" << ba << "\t" << bf << "\t"
				<< rt*rt << "\t" << ra*ra << "\t" << rf*rf << "\t"
				<< sdbt  << "\t" << sdba << "\t" << sdbf << endl;
	}
	return(1);
}



// 	Q: Valores de Q posibles 
//	pixval: Measure
//  fileout: Archivo de Output
//
//  Warning: assume dimX==dimY and minBoxSize=4
//
int MultifractalSBA(simplmat <double> &pixval,simplmat <double> &q, char * outFile, 
	int minBoxSize, int maxBoxSize, int numBoxSizes, char normalize, char * ident)
{
//	simplmat <double> piQ;
//	simplmat <double> piHat;
	simplmat <double> box;
	simplmat <double> tauQ;
	simplmat <double> alphaQ;
	simplmat <double> fQ;
	simplmat <int> boxIni;
	simplmat <int> boxFin;
	
	double sumAlphaQ, sumFQ,cnt,qT,piQT,piHatT,tauQT;

	int  boxSize, iq, i, iRow, iCol, ix, iy, qNum, xDim, yDim;
	int  actBoxSize=0,yResto=0,xResto=0;

	qNum = q.getRows();
	xDim = pixval.getRows();
	yDim = pixval.getCols();

	cnt=0;
   	double tot=0;

	for(iy=0; iy<yDim; iy++)
		for(ix=0; ix<xDim; ix++)
			cnt+=pixval(ix,iy);

    if(toupper(normalize)=='S')
    {
		for(iy=0; iy<yDim; iy++)
			for(ix=0; ix<xDim; ix++)
			{
				pixval(ix,iy)/=cnt;
				tot+=pixval(ix,iy);
			}
	}		
	else if(toupper(normalize)=='D')
	{
		cnt+=xDim*yDim;
		for(iy=0; iy<yDim; iy++)
			for(ix=0; ix<xDim; ix++)
			{
				pixval(ix,iy)=(pixval(ix,iy)+1)/cnt;					
				tot+=pixval(ix,iy);
			}
	}


	if( maxBoxSize > yDim/2 || maxBoxSize > xDim/2  )
		maxBoxSize = (xDim<yDim ? xDim : yDim)/2;

//
//	Calculo de boxes con intervalo constante en escala logaritmica
//
/*
	double deltaBoxSize = (log10(static_cast<double>(maxBoxSize)) - log10(static_cast<double>(minBoxSize)))/ numBoxSizes;
	box.resize(numBoxSizes);

	box(0)= minBoxSize;
	for(i=1; i<numBoxSizes; i++)
	{
		box(i)=pow10(log10(box(i-1))+deltaBoxSize);
		if(box(i) > maxBoxSize )
		{
			numBoxSizes = i;
			break;
		}
		
	}
	for(i=0; i<numBoxSizes; i++)
		box(i) = ceil(box(i));
		
	for(i=numBoxSizes-1; i>1; i--)
	{
		if( box(i) == box(i-1) )
		{
			for(int ii=i; ii<numBoxSizes-1; ii++)
				box(ii) = box(ii+1);
			numBoxSizes--;
		}
	}
*/
	// Calculo de boxes como potencias de 2
	//
	box.resize(numBoxSizes);
   	int finDimX=0,finDimY=0;

	for(finDimX=0,i=0;finDimX<minBoxSize;finDimX=pow(2.0,i))
    	i++;
    int minPot=i;
    
	box(0)= pow(2.0,minPot);
	for(i=1; i<numBoxSizes; i++)
	{
		box(i)=pow(2.0,minPot+i);
		if( box(i)>maxBoxSize )
		{
			numBoxSizes = i;
			break;
		}
	}

	for(finDimX=0;finDimX<xDim;finDimX=pow(2.0,i))
    	i++;
	finDimX=pow(2.0,i-1);
	for(finDimY=0,i=1;finDimY<yDim;finDimY=pow(2.0,i))
    	i++;
	finDimY=pow(2.0,i-1);

	int numRep;

	alphaQ.resize(numBoxSizes,qNum,0.0);
	tauQ.resize(numBoxSizes,qNum,0.0);
	fQ.resize(numBoxSizes,qNum,0.0);
	boxIni.resize(4,2,0);
	boxFin.resize(4,2,0);
	
	for(boxSize=0; boxSize<numBoxSizes; boxSize++)
	{
		for(iq=0;iq<qNum; iq++)
		{
			qT = q(iq);

			actBoxSize = box(boxSize);
			yResto = yDim - finDimY;
			xResto = xDim - finDimX;

			boxIni.fill(0);
			boxFin.fill(finDimX);
			boxFin(0,1)=boxFin(1,1)=boxFin(2,1)=boxFin(3,1)=finDimY;
			numRep=1;

			if( xResto>0 && yResto>0  )
			{
				numRep=4;
				boxIni(1,0)= xResto;
				boxIni(3,0)= xResto;
				boxIni(2,1)= yResto;
				boxIni(3,1)= yResto;
				boxFin(1,0)= finDimX+xResto;
				boxFin(3,0)= finDimX+xResto;
				boxFin(2,1)= finDimY+yResto;
				boxFin(3,1)= finDimY+yResto;

			}
			if( xResto>0  && yResto==0)
			{
				numRep=2;
				boxIni(1,0)= xResto;
				boxFin(1,0)= finDimX+xResto;
			}
			if( yResto>0  && xResto==0)
			{
				numRep=2;
				boxIni(1,1)= yResto;
				boxFin(1,1)= finDimY+yResto;
			}

			for(int rep=0;rep<numRep;rep++)
			{
				piQT=0;
				tauQT=0;
				
				for(iRow=boxIni(rep,1); iRow <= boxFin(rep,1)-actBoxSize; iRow+=actBoxSize )
					for(iCol=boxIni(rep,0); iCol <= boxFin(rep,0)-actBoxSize; iCol+=actBoxSize )
					{
						cnt = 0.0;

//						window examination
						for(iy=iRow; iy<iRow+actBoxSize; iy++)
							for(ix=iCol; ix<iCol+actBoxSize; ix++)
								cnt+=pixval(ix,iy);
						if( cnt>0.0 )
							piQT+=pow(cnt,qT);
						
					}
				sumAlphaQ=0;
				sumFQ=0;
				piHatT=0;

				if( piQT > 0.0 )
					tauQT=log10(piQT);
				else
				{
					goto AFTERSUM;
				}
		
				// To do AlphaQ and FQ

				//	window movement
				for(iRow=boxIni(rep,1); iRow <= boxFin(rep,1)-actBoxSize; iRow +=actBoxSize )
					for(iCol=boxIni(rep,0); iCol <=boxFin(rep,0)-actBoxSize; iCol +=actBoxSize )
					{

						cnt=0.0;

						//	window examination
						for(iy=iRow; iy<iRow+actBoxSize; iy++)
							for(ix=iCol; ix<iCol+actBoxSize; ix++)
								cnt+=pixval(ix,iy);

						if( cnt>0.0 )
						{
							piHatT=pow(cnt,qT)/piQT;
							sumAlphaQ+=piHatT*log10(cnt);
							sumFQ+=piHatT*log10(piHatT);
						}
						
//						if( (piQT != 0.0) && (cnt!=0.0) )
//							piHatT=pow(cnt,qT)/piQT;
//						else
//							piHatT=0;
//
//						if( cnt > 0.0 )
//							sumAlphaQ+=piHatT*log10(cnt);
//						
//						if(piHatT > 0.0)
//							sumFQ+=piHatT*log10(piHatT);
					}

				AFTERSUM:
				alphaQ(boxSize,iq)+=sumAlphaQ;
				fQ(boxSize,iq)+=sumFQ;
				tauQ(boxSize,iq)+=tauQT;
			}

			tauQ(boxSize,iq)/=static_cast<double>(numRep);
			alphaQ(boxSize,iq)/=static_cast<double>(numRep);
			fQ(boxSize,iq)/=static_cast<double>(numRep);
			
		}
	}

	ofstream fb;
	int privez=0;

	fb.open( outFile, ios::in );
	if(!fb )
		privez=1;
	fb.close();
	fb.open( outFile, ios::app );
	if( !fb )
		{
		cerr << "Cannot open patch stats file: " << outFile << endl;
		return 0;
		}

	if( privez )
		fb << "File\tq\tTau\talfa\tf(alfa)\tR-Tau\tR-alfa\tR-f\tSD-Tau\tSD-alfa\tSD-f" << endl;

	// Regression of the log - log variables
	//
	//
	for(i=0; i<qNum; i++)
	{
		double sumx=0,
			sumyt = 0,
			sumya = 0,
			sumyf = 0,
			sumxyt = 0,
			sumxya = 0,
			sumxyf = 0,
			sumxsq = 0,
			sumysqt = 0,
			sumysqa = 0,
			sumysqf = 0,
			sdbt=0,
			sdba=0,
			sdbf=0;

		for(boxSize=0; boxSize<numBoxSizes; boxSize++)
		{
			double x = log10(box(boxSize));
			double yt = tauQ(boxSize,i);
			double ya = alphaQ(boxSize,i);
			double yf = fQ(boxSize,i);
			
			sumx += x;
			sumyt += yt;
			sumya += ya;
			sumyf += yf;
			sumxyt += x * yt;
			sumxya += x * ya;
			sumxyf += x * yf;
			sumxsq += x * x;
			sumysqt += yt * yt;
			sumysqa += ya * ya;
			sumysqf += yf * yf;
		}
		double xbar  = sumx/numBoxSizes ;
		double ybart = sumyt/numBoxSizes ;
		double ybara = sumya/numBoxSizes ;
		double ybarf = sumyf/numBoxSizes ;
	
		double bt = (sumxyt - numBoxSizes *xbar*ybart)/(sumxsq - numBoxSizes *xbar*xbar);
		double ba = (sumxya - numBoxSizes *xbar*ybara)/(sumxsq - numBoxSizes *xbar*xbar);
		double bf = (sumxyf - numBoxSizes *xbar*ybarf)/(sumxsq - numBoxSizes *xbar*xbar);
	//	a = ybar - (b*xbar);
		double rt = (sumxyt - numBoxSizes *xbar*ybart)/sqrt((sumxsq - numBoxSizes  * xbar * xbar) * (sumysqt - numBoxSizes *ybart*ybart));
		double ra = (sumxya - numBoxSizes *xbar*ybara)/sqrt((sumxsq - numBoxSizes  * xbar * xbar) * (sumysqa - numBoxSizes *ybara*ybara));
		double rf = (sumxyf - numBoxSizes *xbar*ybarf)/sqrt((sumxsq - numBoxSizes  * xbar * xbar) * (sumysqf - numBoxSizes *ybarf*ybarf));
		if( numBoxSizes  > 4 )
		{
			sdbt = sqrt((1-rt*rt)*(sumysqt-numBoxSizes*ybart*ybart))/((numBoxSizes -4)*(sumxsq-numBoxSizes *xbar*xbar));
			sdba = sqrt((1-ra*ra)*(sumysqa-numBoxSizes*ybara*ybara))/((numBoxSizes -4)*(sumxsq-numBoxSizes *xbar*xbar));
			sdbf = sqrt((1-rf*rf)*(sumysqf-numBoxSizes*ybarf*ybarf))/((numBoxSizes -4)*(sumxsq-numBoxSizes *xbar*xbar));
		}
		fb << ident << "\t" << q(i) << "\t" << bt << "\t" << ba << "\t" << bf << "\t"
				<< rt*rt << "\t" << ra*ra << "\t" << rf*rf << "\t"
				<< sdbt  << "\t" << sdba << "\t" << sdbf << endl;
	}
	return(1);
}
