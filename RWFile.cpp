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
#include "RWFile.h"
#include <fstream>
#include <iomanip>
#include <string>
//#include "fortify.h"

using namespace std;

int RWFile::WriteIdrisi( const char * fname, simplmat<float>& data)
{
	//int i,j;

	string dname,iname;
	iname = fname;
    string::size_type pos=0;

	if( (pos=iname.find(".rst")) == string::npos )
    {
		iname += ".rst";
		dname = fname;
		dname += ".rdc";
    }
	else
    {
    	dname = iname.substr(0,pos) + ".rdc";
    }
	

	ofstream sav( iname.c_str(), ios::binary | ios::out );
	if(!sav)
	{
		cerr << "Cannot open img file: " << iname.c_str() << endl;
		return 0;
	}
	DimY = data.getCols();
	DimX = data.getRows();
	
	float * f = data.pointer() ;
	
	sav.write((char *)f,DimX*DimY*4);

    sav.close();
    
	//
    // Encuentra Maximo y minimo
    //
	int dy,dx;
	float max=0,min=0;
	for(dy=0;dy<DimY; dy++)
		for(dx=0;dx<DimX; dx++)
		{
			if( data(dx,dy) > max )
				max = data(dx,dy);
			else
				if(data(dx,dy) < min )
					min = data(dx,dy);
		}
		
	sav.open( dname.c_str(), ios::out );
	if(!sav)
	{
		cerr << "Cannot open doc file: " << dname.c_str() << endl;
		return 0;
	}

	sav << "file format : IDRISI Raster A.1" << endl;
	sav << "file title  : "  <<  iname.c_str() << " BI File " << endl;
	sav << "data type   : real" << endl;
	sav << "file type   : binary" << endl;
	sav << "columns     : " << DimX << endl;
	sav << "rows        : " << DimY << endl;
	sav << "ref. system : plane" << endl;
	sav << "ref. units  : m" << endl;
	sav << "unit dist.  : 1.0000000" << endl;
	sav << "min. X      : 1" << endl;
	sav << "max. X      : " << DimX << endl;
	sav << "min. Y      : 1" << endl;
	sav << "max. Y      : " << DimY << endl;
	sav << "pos'n error : unknown" << endl;
	sav << "resolution  : unknown" << endl;
	sav << "min. value  : "<< min << endl;
	sav << "max. value  : "<< max << endl;
	sav << "display min : "<< min << endl;
	sav << "display max : "<< max << endl;
	sav << "value units : unspecified" << endl;
	sav << "value error : unknown" << endl;
	sav << "flag value  : none" << endl;
	sav << "flag def'n  : none" << endl;
	sav << "legend cats : 0" << endl;

   sav.close();


/*	for(i=0; i<DimY; i++)
	{
		for(j=0;j<DimX;j++)
		{
			sav<< data(j,i) << "\t";
		}
		sav << endl;
	}
	sav << endl;
*/
	return 1;

}

int RWFile::WriteIdrisi( const char * fname, simplmat<int>& data)
{
	//int i,j;

	string dname,iname;
	iname = fname;

    string::size_type pos=0;

	if( (pos=iname.find(".rst")) == string::npos )
    {
		iname += ".rst";
		dname = fname;
		dname += ".rdc";
    }
	else
    {
    	dname = iname.substr(0,pos) + ".rdc";
    }
	




	ofstream sav( iname.c_str(), ios::binary | ios::out );
	if(!sav)
	{
		cerr << "Cannot open img file: " << iname.c_str() << endl;
		return 0;
	}
	DimY = data.getCols();
	DimX = data.getRows();
	
//	int * f = data.pointer() ;
	
	int dy,dx,eint=0;
	for(dy=0;dy<DimY; dy++)
		for(dx=0;dx<DimX; dx++)
		{
			eint = data(dx,dy);
			sav.write(reinterpret_cast<char *>(&eint),2);
			}
//	sav.write((char *)f,DimX*DimY*2);

    sav.close();
    
	//
    // Encuentra Maximo y minimo
    //
	int max=0,min=0;
	for(dy=0;dy<DimY; dy++)
		for(dx=0;dx<DimX; dx++)
		{
			if( data(dx,dy) > max )
				max = data(dx,dy);
			else
				if(data(dx,dy) < min )
					min = data(dx,dy);
		}
		
	sav.open( dname.c_str(), ios::out );
	if(!sav)
	{
		cerr << "Cannot open doc file: " << dname.c_str() << endl;
		return 0;
	}

   sav << "file format : IDRISI Raster A.1" << endl;
   sav << "file title  : "  <<  iname.c_str() << " BI File " << endl;
	sav << "data type   : integer" << endl;
	sav << "file type   : binary" << endl;
	sav << "columns     : " << DimX << endl;
	sav << "rows        : " << DimY << endl;
	sav << "ref. system : plane" << endl;
	sav << "ref. units  : m" << endl;
	sav << "unit dist.  : 1.0000000" << endl;
	sav << "min. X      : 1" << endl;
	sav << "max. X      : " << DimX << endl;
   sav << "min. Y      : 1" << endl;
	sav << "max. Y      : " << DimY << endl;
	sav << "pos'n error : unknown" << endl;
	sav << "resolution  : unknown" << endl;
	sav << "min. value  : "<< min << endl;
	sav << "max. value  : "<< max << endl;
	sav << "display min : "<< min << endl;
	sav << "display max : "<< max << endl;
	sav << "value units : unspecified" << endl;
	sav << "value error : unknown" << endl;
	sav << "flag value  : none" << endl;
	sav << "flag def'n  : none" << endl;
	sav << "legend cats : 0" << endl;

   sav.close();


	return 1;
}



	
