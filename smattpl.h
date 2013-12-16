#ifndef __SMATTPL_HPP
#define __SMATTPL_HPP

#include <iostream>
#include <cstdlib>

#ifdef __GNUC__
//#include "fortify.h"
#endif

#define INTEGER int
using namespace std;

template <class Type> class simplmat
	{
	INTEGER rows,cols;
	Type *x;
	public:
	simplmat() {rows=0; cols=0; x=NULL; };
	~simplmat() { if( x!=NULL) delete []x; };

	simplmat(INTEGER r, INTEGER c=1);

	INTEGER getRows() {return rows;};
	INTEGER getCols() {return cols;};
	
	void resize(INTEGER r, INTEGER c=1);
	void resize(INTEGER r, INTEGER c, Type fval);
	Type & elem(INTEGER i, INTEGER j=0)
		{
	#ifdef RANGE_CHECKING
		if ( i < 0	||	rows <= i  ||  j < 0	||	cols <= j )
			error("Error, simplmat, check2, Index out of range");
	#endif
		return x[i + rows*j ];
//		return x[i*cols + j ];
		};

	Type & operator()(INTEGER i){ return elem(i,0); };

	Type & operator()(INTEGER i, INTEGER j){ return elem(i,j); };

	void fill(Type fval);
	void error(const char* msg);
	Type * pointer() { return x; };

	simplmat(simplmat<Type> const& c);
	
	};

template <class Type> simplmat<Type>::simplmat(simplmat<Type> const& c)
{
	rows = c.rows;
	cols = c.cols;
	x = new Type[rows*cols];
	if (x == 0) error("Error, simplmat, simplmat, Operator new failed");
	memcpy(x,c.x,rows*cols*sizeof(Type) );
}


template <class Type> simplmat<Type>::simplmat(INTEGER r, INTEGER c)
	{
	if (r<=0) error("Error, simplmat, simplmat, Number of rows not positive");
	if (c<=0) error("Error, simplmat, simplmat, Number of columns not positive");
	rows=r;
	cols=c;
	x = new Type[r*c];
	if (x == 0) error("Error, simplmat, simplmat, Operator new failed");
	}

template <class Type> void simplmat<Type>::resize(INTEGER r, INTEGER c)
	{
	if( r==rows && c==cols )
		return;

	if (x != NULL)
		delete[] x;
	rows=r;
	cols=c;
	x = new Type[r*c];
	if (x == NULL)
		error ("Error, simplmat, resize, Operator new failed");
	}

template <class Type> void simplmat<Type>::resize(INTEGER r, INTEGER c, Type fval)
	{
	resize(r,c);
	fill(fval);
	}

template <class Type> void simplmat<Type>::error(const char* msg)
	{
	cerr << msg << "\n";;
	exit(1);
	}

template <class Type> void simplmat<Type>::fill(Type fval)
	{
	INTEGER len = rows*cols;
	Type* top = &(x[len]);
	Type* t = x;
	while (t < top) *t++ = fval;
	}

#endif
