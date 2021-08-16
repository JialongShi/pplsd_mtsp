#pragma once

#include "problem.h"

#include <float.h> 

using namespace std;

#define ILLEGAL_COS 2

class solution
{
public:
	vector<int> x;
	vector<double> fitness;

	bool explored;
	double scalarFunctionValue;
	double cosAngleToWV;
	bool inBoundFlag;

public:
	solution();

	void randomInitialize(const problem& inst);

	bool calculateFitness(const problem& inst);

	bool judgeDominate(const solution &sol2) const;

	bool judgeSameSol(const solution &sol2) const;

	bool judgeSameFit(const vector<double> &fit2) const;
};




