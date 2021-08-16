#include "solution.h"

using namespace std;

solution::solution() {
	explored = false;
	scalarFunctionValue = -DBL_MAX;
	cosAngleToWV = ILLEGAL_COS;
	inBoundFlag = false;
}

void solution::randomInitialize(const problem& inst) {
	
	int n = inst.n;
	x.resize(n);
	for (int i = 0; i < n; i++) {
		x[i] = i;
	}
	int temp;
	int randPos;
	for (int i = 0; i < n; ++i) {
		randPos = rand() % (n - i) + i;
		temp = x[i];
		x[i] = x[randPos];
		x[randPos] = temp;
	}
	calculateFitness(inst);
}


bool solution::calculateFitness(const problem& inst) {
	
	int n = inst.n;
	int m = inst.m;
	if (x.size() != n) {
		return false;
	}
	fitness.resize(m);

	int P, Q;
	for (int mi = 0; mi < m; mi++) {
		fitness[mi] = 0;
		for (int i = 0; i < n; i++) {
			P = x[i];
			Q = x[(i == (n - 1)) ? 0 : (i + 1)];
			fitness[mi] += inst.dist[mi][P][Q];
		}
	}
	return true;
}



bool solution::judgeDominate(const solution &sol2) const {
	//verify whether fitness vector A dominate fitness vector B
	if (fitness.size() != sol2.fitness.size()) {
		return false;
	}
	bool allFlag, anyFlag;

	//if A dominate B, means all(A[i]>=B[i]) && any(A[i]>B[i])

	allFlag = true;
	anyFlag = false;

	for (int mi = 0; mi < fitness.size(); mi++) {
		if (fitness[mi] > sol2.fitness[mi]) {
			anyFlag = true;
		}
		if (fitness[mi] < sol2.fitness[mi]) {
			allFlag = false;
			break;
		}
	}
	return allFlag && anyFlag;
}

bool solution::judgeSameSol(const solution &sol2) const {
	if (fitness.size() != sol2.fitness.size() || x.size() != sol2.x.size()) {
		return false;
	}
	int m = fitness.size();
	int n = x.size();

	for (int mi = 0; mi < m; mi++) {
		if (fitness[mi] != sol2.fitness[mi]) {
			return false;
		}
	}
	bool sameFlag = true;

	int posInSol2 = 0;
	while (x[0] != sol2.x[posInSol2]) {
		posInSol2++;
	}

	if (x[1] == sol2.x[(posInSol2 == (n - 1)) ? 0 : (posInSol2 + 1)]) {//same direction
		for (int i = 2; i < n; i++) {
			if (x[i] != sol2.x[(posInSol2 == (n - 1)) ? 0 : (posInSol2 + i)]) {
				sameFlag = false;
				break;
			}
		}

	}
	else if (x[n - 1] == sol2.x[(posInSol2 == (n - 1)) ? 0 : (posInSol2 + 1)]) {//different direction
		for (int i = 2; i < n; i++) {
			if (x[n - i] != sol2.x[(posInSol2 == (n - 1)) ? 0 : (posInSol2 + i)]) {
				sameFlag = false;
				break;
			}
		}
	}
	else {
		sameFlag = false;
	}

	return sameFlag;
}

bool solution::judgeSameFit(const vector<double> &fit2) const {
	if (fitness.size() != fit2.size()) {
		return false;
	}
	int m = fitness.size();

	for (int mi = 0; mi < m; mi++) {
		if (fitness[mi] != fit2[mi]) {
			return false;
		}
	}

	return true;
}








