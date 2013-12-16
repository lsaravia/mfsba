#include "mf.h"
#include "RWFile.h"

int main(int argc, char * argv[])
{
	RWFile file;
	simplmat <double> data;	
	simplmat <double> q;	

	if( argc < 6)
	{
		cerr << "Multifractal estimation of multispecies 2D distributions" << endl;
		cerr << "Takes the most abundant species and assign it the number 1" << endl;
		cerr << "the following abundant species and assign it the number 2" << endl;
		cerr << "and the same for all remaining species" << endl;
		cerr << "The Dq,F(alpha),Alpha estimation use box counting and the canonical method" << endl;
		cerr << "Chhabra, Meneveau, Jensen & Sreenivasan, Phys.Rev.A 40,5284 (1989)" << endl << endl;

		cerr << "Usage: mf inputFile qFile minBox maxBox numBoxSizes option" << endl
        	 << "       option N: Not normalize measure" << endl
        	 << "              S: Normalize measure" << endl
        	 << "              D: ADD 1 and normalize measure" << endl
        	 << "              A: Save the reordered distribution" << endl; 
        	 
        
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
    
    char norm=argv[6][0];
    	
	simplmat <double> newdata(data);	

	if(MultispeciesReordering(data,newdata))
		MultifractalSBA(newdata, q,argv[1] ,minBox, maxBox, deltaBox, norm);

    if(argv[6][0]=='A')
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

