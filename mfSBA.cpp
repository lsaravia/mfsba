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
#include <vector>
#include <unordered_set>
#include <unordered_map>

using namespace std;

// 	q: q exponent for the multifractal estimation 
//	pixval: Measure
//  outFile: Archivo de Output
//
int MultifractalSBA(simplmat <double> &pixval,simplmat <double> &q, char * outFile,
	int minBoxSize, int maxBoxSize, int numBoxSizes, char normalize)
{
	simplmat <double> box;
	simplmat <double> tauQ;
	simplmat <double> alphaQ;
	simplmat <double> fQ;
	
	if(normalize=='E')
		standardBoxCountSAD(pixval,q, minBoxSize, maxBoxSize, numBoxSizes, normalize,
    									box,tauQ, alphaQ, fQ);
	else
		standardBoxCount(pixval,q, minBoxSize, maxBoxSize, numBoxSizes, normalize,
    									box,tauQ, alphaQ, fQ, &winMovSum);

	int qNum = q.getRows();
	int i,boxSize;

	string tFileName("t.");
	tFileName+=outFile;
	ofstream tFile(tFileName.c_str());
   
	string aFileName("a.");
	aFileName+= outFile;
	ofstream aFile(aFileName.c_str());
   
	string fFileName("f.");
	fFileName+= outFile;
	ofstream fFile(fFileName.c_str());

	tFile << "BoxSize" << "\t" << "LogBox";
	aFile << "BoxSize" << "\t" << "LogBox";
	fFile << "BoxSize" << "\t" << "LogBox";
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
		// Use areas instead of side for SAD multifractals
		if(normalize=='E')
			box(boxSize) = pow(box(boxSize),2);

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

	simplmat <outRegress> outR(qNum);

	loglogRegress(q,numBoxSizes,box,tauQ,alphaQ,fQ,outR);

	for(i=0; i<qNum; i++)
	{
		sFile << q(i) << "\t" << outR(i).bt << "\t" << outR(i).ba << "\t" << outR(i).bf << "\t"
				<< outR(i).rt2 << "\t" << outR(i).ra2 << "\t" << outR(i).rf2 << "\t"
				<< outR(i).sdbt  << "\t" << outR(i).sdba << "\t" << outR(i).sdbf << endl;
	}
	return(1);
}


//  Multifractal estimation with a one-line output used for calling from other programs
//
// 	q: q exponent for the multifractal estimation 
//	pixval: Measure
//  ident: identification for the output 
//
//
int MultifractalSBA(simplmat <double> &pixval,simplmat <double> &q, char * outFile, 
	int minBoxSize, int maxBoxSize, int numBoxSizes, char normalize, char * ident)
{
	simplmat <double> box;
	simplmat <double> tauQ;
	simplmat <double> alphaQ;
	simplmat <double> fQ;

	if(normalize=='E')
		standardBoxCountSAD(pixval,q, minBoxSize, maxBoxSize, numBoxSizes, normalize,
    									box,tauQ, alphaQ, fQ);
	else
		standardBoxCount(pixval,q, minBoxSize, maxBoxSize, numBoxSizes, normalize,
    									box,tauQ, alphaQ, fQ, &winMovSum);

	int qNum = q.getRows();
	ofstream fb;
	int privez=0,i;

	fb.open( outFile, ios::in );
	if(!fb )
		privez=1;
	fb.close();
	fb.open( outFile, ios::app );
	if( !fb )
		{
		cerr << "Cannot open multifractal output file: " << outFile << endl;
		return 0;
		}

	if( privez )
		fb << "File\tq\tTau\talfa\tf(alfa)\tR-Tau\tR-alfa\tR-f\tSD-Tau\tSD-alfa\tSD-f" << endl;

	simplmat <outRegress> outR(qNum);

	// Use areas instead of side for SAD multifractals
	if(normalize=='E')
	{
		for(i=0; i<numBoxSizes; i++)
			box(i) = pow(box(i),2);
	}

	loglogRegress(q,numBoxSizes,box,tauQ,alphaQ,fQ,outR);

	for(i=0; i<qNum; i++)
	{
		fb << ident << "\t" << q(i) << "\t" << outR(i).bt << "\t" << outR(i).ba << "\t" << outR(i).bf << "\t"
				<< outR(i).rt2 << "\t" << outR(i).ra2 << "\t" << outR(i).rf2 << "\t"
				<< outR(i).sdbt  << "\t" << outR(i).sdba << "\t" << outR(i).sdbf << endl;
	}
	return(1);
}


int loglogRegress(simplmat <double> &q,int &numBoxSizes, simplmat <double> &box,
	simplmat <double> &tauQ, simplmat <double> &alphaQ, simplmat <double> &fQ, simplmat <outRegress> &oR)
{
	int qNum = q.getRows();

	// Regression of the log - log variables
	//
	//

	for(int i=0; i<qNum; i++)
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

		for(int boxSize=0; boxSize<numBoxSizes; boxSize++)
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

		oR(i).bt = bt;
		oR(i).ba = ba;
		oR(i).bf = bf;
		oR(i).rt2 = rt*rt;
		oR(i).ra2 = ra*ra;
		oR(i).rf2 = rf*rf;
		oR(i).sdbt = sdbt;
		oR(i).sdba = sdba;
		oR(i).sdbf = sdbf;
	}
	return(1);
}

double winMovSum(simplmat <double> &pixval,const int &rowIni,const int &rowEnd,const int &colIni,const int &colEnd)
{
	double cnt=0;
	for(int iy=rowIni; iy<rowEnd; iy++)
		for(int ix=colIni; ix<colEnd; ix++)
			cnt+=pixval(ix,iy);
	return(cnt);
}

double winMovNumSp(simplmat <double> &px,const int &rowIni,const int &rowEnd,const int &colIni,const int &colEnd)
{
   	unordered_set<int> distinct_container;

	for(int iy=rowIni; iy<rowEnd; iy++)
		for(int ix=colIni; ix<colEnd; ix++)
	   	{
   			int curr_int = static_cast <int>(px(ix,iy));
	    	distinct_container.insert(curr_int);
   		}

   	return distinct_container.size();
} 

int standardBoxCount(simplmat <double> &pixval,simplmat <double> &q, int &minBoxSize, 
    int &maxBoxSize, int &numBoxSizes, char &normalize,
    simplmat <double> &box, simplmat <double> &tauQ, simplmat <double> &alphaQ, simplmat <double> &fQ,
    double (*winMov)(simplmat <double> &pixval,const int &rowIni,const int &rowEnd,const int &colIni,const int &colEnd))
{
	double sumAlphaQ, sumFQ,cnt,qT,piQT,piHatT,tauQT;

	int  boxSize, iq, i, iRow, iCol, ix, iy, qNum, xDim, yDim;
	int  yResto=0,xResto=0,actBoxSize=0;

	simplmat <int> boxIni;
	simplmat <int> boxFin;

	qNum = q.getRows();
	xDim = pixval.getRows();
	yDim = pixval.getCols();

	cnt=0.0;
   	normalize = toupper(normalize);

   	if(normalize=='S') 
   	{
		for(iy=0; iy<yDim; iy++)
			for(ix=0; ix<xDim; ix++)
				cnt+=pixval(ix,iy);
	
		for(iy=0; iy<yDim; iy++)
			for(ix=0; ix<xDim; ix++)
				pixval(ix,iy)/=cnt;
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

	for(finDimX=0;finDimX<=xDim;finDimX=pow(2.0,i))
    	i++;
	finDimX=pow(2.0,i-1);
	for(finDimY=0,i=1;finDimY<=yDim;finDimY=pow(2.0,i))
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
				piQT=0.0;
				tauQT=0.0;
				sumAlphaQ=0.0;
				sumFQ=0.0;
				piHatT=0.0;
				
				for(iRow=boxIni(rep,1); iRow <= boxFin(rep,1)-actBoxSize; iRow+=actBoxSize )
					for(iCol=boxIni(rep,0); iCol <= boxFin(rep,0)-actBoxSize; iCol+=actBoxSize )
					{
//						window examination
						cnt = winMov(pixval,iRow,iRow+actBoxSize,iCol,iCol+actBoxSize);
						if( cnt>0.0 )
							piQT+=pow(cnt,qT);
						
					}
				

				if( piQT > 0.0 )
				{
					tauQT=log10(piQT);
					// To do AlphaQ and FQ

					//	window movement
					for(iRow=boxIni(rep,1); iRow <= boxFin(rep,1)-actBoxSize; iRow +=actBoxSize )
						for(iCol=boxIni(rep,0); iCol <=boxFin(rep,0)-actBoxSize; iCol +=actBoxSize )
						{

							cnt = winMov(pixval,iRow,iRow+actBoxSize,iCol,iCol+actBoxSize);

							if( cnt>0.0 )
							{
								piHatT=pow(cnt,qT)/piQT;
								sumAlphaQ+=piHatT*log10(cnt);
								sumFQ+=piHatT*log10(piHatT);
							}
							
						}
				}

				alphaQ(boxSize,iq)+=sumAlphaQ;
				fQ(boxSize,iq)+=sumFQ;
				tauQ(boxSize,iq)+=tauQT;
			}
			
			tauQ(boxSize,iq)/=static_cast<double>(numRep);
			alphaQ(boxSize,iq)/=static_cast<double>(numRep);
			fQ(boxSize,iq)/=static_cast<double>(numRep);
			
		}
	}
	return(1);
}

int standardBoxCountSAD(simplmat <double> &pixval,simplmat <double> &q, int &minBoxSize, 
    int &maxBoxSize, int &numBoxSizes, char &normalize,
    simplmat <double> &box, simplmat <double> &tauQ, simplmat <double> &alphaQ, simplmat <double> &fQ)
{
	double sumAlphaQ, sumFQ,cnt,qT,piQT,tauQT;

	int  boxSize, iq, i, iRow, iCol, ix, iy, qNum, xDim, yDim;
	int  yResto=0,xResto=0,actBoxSize=0;

	simplmat <int> boxIni;
	simplmat <int> boxFin;

	qNum = q.getRows();
	xDim = pixval.getRows();
	yDim = pixval.getCols();

	cnt=0.0;
   	normalize = toupper(normalize);

   	if(normalize=='S') 
   	{
		for(iy=0; iy<yDim; iy++)
			for(ix=0; ix<xDim; ix++)
				cnt+=pixval(ix,iy);
	
		for(iy=0; iy<yDim; iy++)
			for(ix=0; ix<xDim; ix++)
				pixval(ix,iy)/=cnt;
	}

	
	if( maxBoxSize > yDim/2 || maxBoxSize > xDim/2  )
		maxBoxSize = (xDim<yDim ? xDim : yDim)/2;

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

	for(finDimX=0;finDimX<=xDim;finDimX=pow(2.0,i))
    	i++;
	finDimX=pow(2.0,i-1);
	for(finDimY=0,i=1;finDimY<=yDim;finDimY=pow(2.0,i))
    	i++;
	finDimY=pow(2.0,i-1);

	int numRep;

	alphaQ.resize(numBoxSizes,qNum,0.0);
	tauQ.resize(numBoxSizes,qNum,0.0);
	fQ.resize(numBoxSizes,qNum,0.0);
	boxIni.resize(4,2,0);
	boxFin.resize(4,2,0);
	for(iq=0;iq<qNum; iq++)
	{
		
		for(boxSize=0; boxSize<numBoxSizes; boxSize++)
		{
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

				qT = q(iq);
				piQT=0.0;
				tauQT=0.0;
				sumAlphaQ=0.0;
				sumFQ=0.0;
				int countBoxes=0;

				for(iRow=boxIni(rep,1); iRow <= boxFin(rep,1)-actBoxSize; iRow+=actBoxSize )
					for(iCol=boxIni(rep,0); iCol <= boxFin(rep,0)-actBoxSize; iCol+=actBoxSize )
					{

						// boxes examination
						piQT += winMovSAD(pixval,iRow,iRow+actBoxSize,iCol,iCol+actBoxSize,qT,
							sumAlphaQ,sumFQ);
						countBoxes++;
					}
				// Calculates average
				piQT /= static_cast<double>(countBoxes);
				sumAlphaQ /= static_cast<double>(countBoxes);
				sumFQ /= static_cast<double>(countBoxes);

				if( piQT > 0.0 )
					tauQT=log10(piQT);

				alphaQ(boxSize,iq)+=sumAlphaQ;
				fQ(boxSize,iq)+=sumFQ;
				tauQ(boxSize,iq)+=tauQT;
			}
			// multipling by -1 because formulae have inverted signs!!!! 
			// compare Saravia with Borda-de-Agua
			//
			tauQ(boxSize,iq)/=static_cast<double>(-numRep);
			alphaQ(boxSize,iq)/=static_cast<double>(-numRep);
			fQ(boxSize,iq)/=static_cast<double>(-numRep);
				
		}
	}
	return(1);
}

// This routine assumes one individual by cell!
//
double winMovSAD(simplmat <double> &px,const int &rowIni,const int &rowEnd,const int &colIni,const int &colEnd,const double &qT,
	double &AlphaQ, double &FQ)
{
	double piQT=0,sumSp=0;
	typedef unordered_map<int,unsigned int> CounterMap;
	CounterMap counts;
	for(int iy=rowIni; iy<rowEnd; iy++)
		for(int ix=colIni; ix<colEnd; ix++)
	   	{
			int curr_int = static_cast <int>(px(ix,iy));
			if(curr_int!=0)
				{
				sumSp++;
				CounterMap::iterator i(counts.find(curr_int));
				if (i != counts.end()){
					i->second++;
			   	} else {
			    	counts[curr_int] = 1;
			   }
		   }
		}

	for(auto iter=counts.begin(); iter!=counts.end(); ++iter)
	{
		piQT += pow((iter->second)/sumSp,qT);
	}


	for(auto iter=counts.begin(); iter!=counts.end(); ++iter)
	{
		double cnt = (iter->second)/sumSp;
		double piHatT=pow(cnt,qT)/piQT;
		AlphaQ+=piHatT*log10(cnt);
		FQ+=piHatT*log10(piHatT);
	}
	
	return(piQT);
} 
