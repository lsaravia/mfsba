#include <cmath>
// Randon generator based on numerical recipes

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

typedef unsigned long long int Ullint;
typedef unsigned int Uint;

struct Ranf1 {
	// Fast randon generator based on numerical recipes, the period is 1.8 10^19 .
	Ullint v;
	Ranf1(Ullint j=0) : v(4101842887655102017LL) {
		if( j==0 ) j=time(NULL);
		v ^= j;
		v = int64();
	}

	inline void init(Ullint j){
		if(j>0){
		v = 4101842887655102017LL ^ j;
		v = int64();
		}
	}

	inline Ullint int64() {
		v ^= v >> 21; v ^= v << 35; v ^= v >> 4;
		return v * 2685821657736338717LL;
	}

	inline double doub() { return 5.42101086242752217E-20 * int64(); }

	inline Uint int32() { return (Uint)int64(); }
};


struct Normaldev : Ranf1 {

	double mu,sig;
	Normaldev(double mmu=0, double ssig=1, Ullint i=0)
	: Ranf1(i), mu(mmu), sig(ssig){}
	//Constructor arguments are mu=median, sig= , and a random sequence seed.
	
	inline void init(double mmu, double ssig, Ullint i){
		mu= mmu;
		sig= ssig;
		Ranf1::init(i);
	}

	inline double dev() {
		//Return a normal deviate.
		double u,v,x,y,q;
		do {
			u = doub();
			v = 1.7156*(doub()-0.5);
			x = u - 0.449871;
			y = std::abs(v) + 0.386595;
			q = SQR(x) + y*(0.19600*y-0.25472*x);
		} while (q > 0.27597 && (q > 0.27846 || SQR(v) > -4.*log(u)*SQR(u)));

		return mu + sig*v/u;
	}
};

struct Lognormaldev : Normaldev
{
	double mu,sig;
	Lognormaldev(double mmu, double ssig, Ullint i) :
	Normaldev(0.0,1.0,i), mu(mmu), sig(ssig){}

	inline double dev()
	{
		return exp(mu+sig*Normaldev::dev());
	}
};