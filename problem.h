#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <string.h>
#include <math.h>

using namespace std;

#define CC_EUCLIDEAN      0
#define CC_EUCLIDEAN_CEIL 1
#define CC_GEOGRAPHIC     2
#define CC_ATT            3

#define GH_PI (3.141592)

class problem {
public:
	string problemName;
	int n;
	int m;
	int	edgeWeightType;  //edge weight type
	vector< vector<double> >	X; //x coord
	vector< vector<double> >	Y; //y coord
	vector< vector< std::vector<int> > > dist; //distance matrix

public:
	//since TSP is a minimization problem, we transfer it to a maximization problem by multipling -1
	bool readProblem_TransferToMaxProblem(const char* filename);
};
