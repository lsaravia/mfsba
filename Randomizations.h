#ifndef RANDOMIZATIONS_H
#define RANDOMIZATIONS_H

#include "smattpl.h"
#include <time.h>
#include <algorithm>
#include "ran.h"
//#include "randlib.h"
//#include "r250.h"

struct randomizePosXY
{
	int x;
	int y;
};

class R250  : Ranf1
{
	public:
	R250() : Ranf1(){
//		int seed = static_cast<long>(time(0));
//		setall(seed,seed+1);
		};
	int operator()(const int &num) 
        {
        	return (int64() % (num));
//          return ignuin(0,num-1);
        }

	double operator()() 
        {
			return doub(); 
//	        return ranf();
        }

};


class Randomizations
{
	R250 rnd;
	public:
	Randomizations() {};

	void Randomize(simplmat<double>& data,simplmat<double> & rdata);
	void Randomize(simplmat<double>& data);
    	
/*	template <class T> void Randomize(T & data){
		T * ptr = data.pointer();
		int max = sizeof(ptr)/sizeof(*ptr);
	    random_shuffle(ptr, ptr + max, rnd);
		};
*/
};
#endif  
