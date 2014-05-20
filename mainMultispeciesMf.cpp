#include "mf.h"
#include "RWFile.h"

int main(int argc, char * argv[])
{
	RWFile file;
	simplmat <double> data;	
	simplmat <double> q;	

	if( argc < 6)
	{
		cerr << "-------------------------------------------------------" << endl 
			<< "Multifractal estimation of multispecies 2D distributions" << endl
			 << "Using species rank surface (SRS) see " << endl 
			 << "http://f1000research.com/articles/3-14/v2" << endl 
			 << "The Dq,F(alpha),Alpha estimation use box counting and the canonical method" << endl
			 << "Chhabra, Meneveau, Jensen & Sreenivasan, Phys.Rev.A 40,5284 (1989)" << endl << endl;

		cerr << "Usage: mf inputFile qFile minBox maxBox numBoxSizes option" << endl
        	 << "       option N: Not normalize measure use SRS" << endl
        	 << "              S: Normalize measure" << endl
        	 << "              A: Save the reordered distribution" << endl 
			 << "              E: Use species number to compute dimensions" << endl	        	 
			 << "                 This option implies N" << endl;	        	 
        
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
   	
   
	if(!file.ReadSeed(argv[2], q))
		exit(1);

	int minBox = atoi(argv[3]);
	int maxBox = atoi(argv[4]);
	int deltaBox = atoi(argv[5]);
    
    string opt=argv[6];
    	
	simplmat <double> newdata(data);	

	if(MultispeciesReordering(data,newdata))
		MultifractalSBA(newdata, q,argv[1] ,minBox, maxBox, deltaBox, opt[0]);

    if(opt[0]=='A')
    {
	    string::size_type pos=0;
		if( (pos=fname.find(".")) == string::npos )
	    {
			fname += "Reord.sed";
	    }
		else
		{
	    	fname = fname.substr(0,pos) + "Reord.sed";
		}

    	if(!file.WriteSeed(fname.c_str(),newdata))
			exit(1);
    }
	return 0;
}

