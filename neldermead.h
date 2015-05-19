#pragma once
#include <vector>


class Vector {
	std::vector<double> _coords;
public:
	Vector(){}
	Vector(std::vector<double> coords) {
		_coords = coords;
	}
	Vector(const Vector& cp) {
		_coords = cp._coords;
	}
	Vector operator*(double a) const {
		auto coords = _coords;
		for(auto it = coords.begin(); it != coords.end(); ++it) {
			*it *= a;
		}
		return Vector(coords);
	}
	Vector operator+(const Vector& arg) const {
		auto coords = arg._coords;
		for(int i = 0; i < _coords.size(); ++i) {
			coords[i] += _coords[i];
		}
		return Vector(coords);
	}
	Vector operator-(const Vector& arg) const {
		auto coords = arg._coords;
		for(int i = 0; i < _coords.size(); ++i) {
			coords[i] -= _coords[i];
		}
		return Vector(coords);
	}
	double get(unsigned i) const {
		return _coords.at(i);
	}

};

class Function
{
protected:
	unsigned _n;			// Number of arguments
public:
	double virtual Calculate(const Vector&) = 0;
	unsigned GetArity() {return _n;}	
};

class NelderMead {
	std::vector<Vector> _simplex;
	std::vector<double> _values;
	Function* _f;
	double _alpha, _beta, _gamma;
	double _arity;
	int minInd;
public:
	NelderMead(Function*, double = 1, double = 0.5, double = 2);
	~NelderMead();
	void step();
	double GetDispersion();
	Vector GetResult(double);		// Gets epsilon
};