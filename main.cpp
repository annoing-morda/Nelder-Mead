#include "neldermead.h"
#include <iostream>
#include <cmath>

class F1 : public Function
{

public:
	F1() {_n = 2;};
	~F1() {};
	double Calculate (const Vector& v) {
		double x = v.get(0) ;
		double y = v.get(1) ;
		return (x)*(x) + (y-2)*(y-2);
	}
};

class F2 : public Function
{

public:
	F2() {_n = 2;};
	~F2() {};
	double Calculate (const Vector& v) {
		double x = v.get(0) ;
		double y = v.get(1) ;
		return -(std::sin(x) + std::sin(y)) * (std::sin(x) + std::sin(y));
	}
};

int main(){
	NelderMead opt (new F1());
	Vector v = opt.GetResult(1e-8);
	return 0;
}