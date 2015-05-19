#include "neldermead.h"
#include <vector>
#include <utility>
#include <cmath>

#include <iostream>

NelderMead::NelderMead(Function* f, double alpha, double beta, double gamma) {
	_f = f;
	_alpha = alpha;
	_beta = beta;
	_gamma = gamma;
	_arity = _f->GetArity();
	minInd = 0;
	_simplex.resize(_arity + 1);
	_values.resize(_arity + 1);
	std::vector<double> coords(_arity);
	for(auto it = coords.begin(); it != coords.end(); ++it) {
		*it = 0;
	}
	_simplex[0] = Vector(coords);
	_values[0] = _f->Calculate(_simplex[0]);
	for(int i = 0; i < _arity; ++i) {
		coords[i] = 1;
		_simplex[i + 1] = Vector(coords);
		_values[i + 1] = _f->Calculate(_simplex[i + 1]);
		coords[i] = 0;
	}
}

NelderMead::~NelderMead() {
	delete _f;
}

void NelderMead::step() {
	unsigned hInd = 0, gInd = 0, lInd = 0;
	for (int i = 1; i < _arity + 1; ++i) {			// Getting xh, xg, xl.
		if (_values[i] > _values[hInd]) {
			gInd = hInd;
			hInd = i;
		} else {
			if (_values[i] > _values[gInd]) {
				gInd = i;
			}
		}
		if (_values[i] < _values[lInd]) {
			lInd = i;
		}
	}
	minInd = lInd;
	Vector xc (_simplex[_arity]);					// Getting center.
	for(int i = 0; i < _arity; ++i) {
		if(i != hInd) {
			xc = xc + _simplex[i];
		}
	}
	xc = xc * (1/(double)_arity);
	Vector xr(xc*(1.0 + _alpha) + _simplex[hInd] *(-_alpha));
	double fr = _f->Calculate(xr);
	if (fr < _values[lInd]) {							// fr < fl
		Vector xe(xc*(1.0 - _gamma) + xr*(_gamma));		// Extension.
		double fe = _f->Calculate(xe);
		if (fe < _values[lInd]) {
			_simplex [hInd] = xe;
			_values[hInd] = fe;
			return;
		}
		_simplex [hInd] = xr;
		_values[hInd] = fr;
		return;
	}
	if (fr < _values[gInd]) {							// fl < fr < fg
		_simplex[hInd] = xr;
		_values[hInd] = fr;
		return;
	}
	if (fr < _values[hInd]) {
		std::swap(xr, _simplex[hInd]);			// ?????
		std::swap(fr, _values[hInd]);			
	}
	Vector xs(_simplex[hInd]*_beta + xc*(1.0 - _beta));	 	// Compression.
	double fs = _f->Calculate(xs);
	if(fs < _values[hInd]) {
		_simplex[hInd] = xs;
		_values[hInd] = fs;
		return;
	}
	for (int i = 0; i < _arity + 1; ++i) {				// Compression of entire simplex
		if (i != lInd) {
			_simplex[i] = (_simplex[i] + _simplex[lInd]) * 0.5;
			_values[i] = _f->Calculate (_simplex [i]);
		}
	}
}

double NelderMead::GetDispersion() {
	double sumR = 0;
	for (int i = 0; i < _arity + 1; ++i) {
		double r = 0;
		for (int j = 0; j < _arity; ++j) {
			r += (_simplex[i].get(j) - _simplex[minInd].get(j))* (_simplex[i].get(j) - _simplex[minInd].get(j));
		}
		sumR += std::sqrt(r);
	}
	return sumR / ((double)_arity + 1);
}

Vector NelderMead::GetResult (double eps) {
	eps = (eps < 0) ? -eps : eps;
	double disp;
	do {
		disp = GetDispersion();
		step();
	} while (disp > eps);
	return _simplex[minInd];
}