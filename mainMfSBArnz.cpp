#include <sstream>
#include <float.h>
#include "mf.h"
#include "RWFile.h"
#include "Randomizations.h"
//#include "Spectral.h"
#include <stdio.h>

int CalcRandomizationMfSBAEnvelope(	const char * inFile, const char * outFile, simplmat <double> &q,const int numSimul, const double probConf=0.05);
void spectMin(double const &pspect, simplmat<double> & pMin, int & l, int nExtr );
void spectMax(double const &pspect, simplmat<double> & pMax, int & l, int nExtr );

int main(int argc, char * argv[])
{
	RWFile file;
	simplmat <double> data;	
	simplmat <double> q;	
	simplmat <double> rdata;
	string nfOut;
	
	if(argc == 5)
	{
		// Lee un archivo de randomization previamente generado 
		// mfSBArnz fileRnz q.sed numsimul P
		//

		nfOut=argv[1];
		nfOut += ".rnz";
		if(!file.ReadSeed(argv[2], q))
			exit(1);
			
		int numSimul = atoi(argv[3]);
		double prob= atof(argv[4]);
		CalcRandomizationMfSBAEnvelope(argv[1],nfOut.c_str(),q,numSimul,prob);
		return 0;
	}
	else if( argc < 8)
	{
		cerr << "Randomize the image and calculates Generalized Dimensions for numSimul times"  << endl;
		cerr << "and calculates a P confidence interval"  << endl;
		cerr << "Multifractal Dq,F(alpha),Alpha estimation using the canonical method" << endl;
		cerr << "Chhabra, Meneveau, Jensen & Sreenivasan, Phys.Rev.A 40,5284 (1989)" << endl << endl;
      	 
        cerr << "Usage: mfSBArnz inputFile outFile qFile minBox maxBox deltaBox option numSimul P" << endl
        	 << "       option N: Not normalize measure" << endl
        	 << "              S: Normalize measure" << endl << endl;

		cerr << "Usage: mfSBArnz fileRnz q.sed numSimul P" << endl;
		exit(1);
	}

   	string fname = argv[1];

	if( fname.find(".rst") != string::npos )
    {
		if(!file.ReadIdrisi(argv[1], data))
			exit(1);
    }
	else if( fname.find(".tif")!=string::npos )
	{
		if(!file.ReadTiff(argv[1], data))
			exit(1);
	}
	else
    {
		if(!file.ReadSeed(argv[1], data))
			exit(1);
    }


	if(!file.ReadSeed(argv[3], q))
		exit(1);

	remove( argv[2] );

	int minBox = atoi(argv[4]);
	int maxBox = atoi(argv[5]);
	int deltaBox = atoi(argv[6]);
	int numSimul = atoi(argv[8]);
	double prob=0.05;
	if (argc>9)
		prob= atof(argv[9]);

    char norm=argv[7][0];
    if(norm=='A') norm='S';


	MultifractalSBA(data, q,argv[2] ,minBox, maxBox, deltaBox, norm,argv[1]);

	Randomizations rz;
	norm='N';			// si los datos están normalizados no los vuelve a normalizar
	int s;
	for(s=0; s<numSimul; s++)
	{
		ostringstream name;
		name << "rnz-" << s << ends;
		rz.Randomize(data);
		//file.WriteSeed(name.str().c_str(),data);
		MultifractalSBA(data, q,argv[2] ,minBox, maxBox, deltaBox, norm,const_cast<char *>(name.str().c_str()));
		name.str("");
	}

	nfOut=argv[2];
	nfOut += ".rnz";
   	
	CalcRandomizationMfSBAEnvelope(argv[2],nfOut.c_str(),q,numSimul,prob);
	return 0;
}

// inputFile outFile qFile minBox maxBox deltaBox option numSimul
int CalcRandomizationMfSBAEnvelope(	const char * inFile, const char * outFile, 
				simplmat <double> &q,const int numSimul, const double probConf)
{
	// read output of randomizations and calculate confidence envelope 
	//
	ifstream in;
	string buff;
	
	in.open(inFile);
	if ( !in )
	{
		cerr << "Cannot open file: " << inFile << endl;
		return false;
	}
	
	int linesBySimul=q.getRows();
	int line=0;
	simplmat <double> Dq(linesBySimul);
	double alfa=0,qq=0,Dqq=0;
	
	int numExtreme = 1+probConf*numSimul/2;
	if( numExtreme <=1 )
	{
		cerr << "numExtreme <= 1" << endl;
		exit(1);
	}
        
	simplmat <double> DqMin(linesBySimul,numExtreme); // Poner el maximo valor posible de
	simplmat <double> DqMax(linesBySimul,numExtreme); // Poner el maximo valor posible de
	DqMax.fill(DBL_MIN);
	DqMin.fill(DBL_MAX);
	
	// linea de titulos
	getline(in,buff);
	
	// Toma los valores reales
	// 
	while( line<linesBySimul )
	{
		getline(in,buff);
		istringstream ins(buff.c_str());
		ins.ignore(256,'\t');
		ins >> qq;
		ins >> Dqq;
		ins >> alfa;
		if(qq!=1)
			Dq(line)=Dqq/(qq-1);
		else
			Dq(line)=alfa;
		
		if( q(line)!=qq )
			exit(0);
		
		line++;
		if(in.eof())
		{
			cerr << "Inconsistent file - " << inFile << endl;
			exit(0);
		}
	}
	int n=0;
	
	while( n<numSimul )
	{
		for(line=0;line<linesBySimul;line++)
		{
			getline(in,buff);
			istringstream ins(buff.c_str());
			ins.ignore(256,'\t');
			ins >> qq;
			ins >> Dqq;
			ins >> alfa;
			if(qq!=1)
				Dqq /=(qq-1);
			else
				Dqq =alfa;
			
			spectMin(Dqq, DqMin,line,numExtreme);
			spectMax(Dqq, DqMax,line,numExtreme);
			if(in.eof())
			{
				cerr << "Inconsistent file - " << inFile << endl;
				exit(0);
			}
		}
		n++;
	}
	ofstream fOut(outFile, ios::out | ios::app );
	fOut << "Randomizations: \t" << numSimul << "\t" << numExtreme << endl;
	fOut << "P:" << "\t" << probConf<< endl;
	fOut << "q\tDq\tDqMax\tDqMin" << endl;
	
	for(line=0;line<linesBySimul;line++)
	{
		int signif = 0;
		fOut << setw(4);
		fOut << q(line) << "\t";
		fOut << Dq(line) << "\t";
		fOut << DqMax(line) << "\t";
		fOut << DqMin(line) << "\t";
		if( Dq(line)<DqMin(line) )
	 			signif = -1;
	 	if( Dq(line)>DqMax(line) )
	 			signif = 1;
	 	fOut << signif << endl;
	}
}

// Mantiene los "nExtr" valores minimos de "pspect"
//
void spectMin(double const &pspect, simplmat<double> & pMin, int & l, int nExtr )
{
	if( pspect < pMin(l,0) )
	{
		for(int k=nExtr-1; k>=0; k--)
		{
			if( pspect < pMin(l,k) )
			{
				for(int kk=1; kk<=k; kk++)
					pMin(l,kk-1) = pMin(l,kk);
				pMin(l,k) = pspect;
				break;
			}
		}
	}
}

void spectMax(double const &pspect, simplmat<double> & pMax, int & l, int nExtr )
{
	if( pspect > pMax(l,0) )
	{
		for(int k=nExtr-1; k>=0; k--)
		{
			if( pspect > pMax(l,k) )
			{
				for(int kk=1; kk<=k; kk++)
					pMax(l,kk-1) = pMax(l,kk);
				pMax(l,k) = pspect;
				break;
			}

		}
	}
}


