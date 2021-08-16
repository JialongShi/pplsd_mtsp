#include "misc.h"

using namespace std;

#if defined(__linux__) || defined(__linux)

void start_timer() {
	ClkTck = sysconf(_SC_CLK_TCK);
	clock_t_start = times(&Old_Time);
}

double return_CPU_time() {
	double cpu_time;
	clock_t_end = times(&New_Time);
	cpu_time = (New_Time.tms_utime - Old_Time.tms_utime) / ((double)ClkTck)
		+ (New_Time.tms_stime - Old_Time.tms_stime) / ((double)ClkTck)
		+ (New_Time.tms_cutime - Old_Time.tms_cutime) / ((double)ClkTck)
		+ (New_Time.tms_cstime - Old_Time.tms_cstime) / ((double)ClkTck);
	return cpu_time;
}

double return_elapsed_time() {
	double real_time;
	clock_t_end = times(&New_Time);
	real_time = (clock_t_end - clock_t_start) / ((double)ClkTck);
	return real_time;
}

#else

void start_timer() {
	time(&Old_Time);
}

double return_elapsed_time() {
	time(&New_Time);
	return (double)difftime(New_Time, Old_Time);
}

#endif

double moveCost2Opt(const int &mi, const int &P, const int &Q, const solution &sol, const problem &inst) {
	int pNext, qNext;
	if (P == Q) {
		return 0;
	}
	int n = inst.n;

	pNext = (P == (n - 1)) ? 0 : (P + 1);
	qNext = (Q == (n - 1)) ? 0 : (Q + 1);

	return inst.dist[mi][sol.x[P]][sol.x[Q]] + inst.dist[mi][sol.x[pNext]][sol.x[qNext]]
		- inst.dist[mi][sol.x[P]][sol.x[pNext]] - inst.dist[mi][sol.x[Q]][sol.x[qNext]];
}

bool move2Opt(const int &P, const int &Q, solution &sol) {
	int n = sol.x.size();
	int cityP = (P <= Q) ? P : Q;
	int cityQ = (P >= Q) ? P : Q;

	int right = (cityP == (n - 1)) ? 0 : (cityP + 1);
	int left = cityQ;
	int orderR, orderL;

	while (right < left) {
		if (right == n)
			right = 0;
		if (left == -1)
			left = n - 1;
		orderR = sol.x[right];
		orderL = sol.x[left];
		sol.x[left--] = orderR;
		sol.x[right++] = orderL;
	}
	return true;
}

double calculateCos_VectorToVector(const vector<double> &V1, const vector<double> &V2) {
	int m = V2.size();
	if (m != V1.size()) return ILLEGAL_COS;

	double V1V2Multiply = 0;
	for (int mi = 0; mi < m; ++mi) {
		V1V2Multiply += V1[mi] * V2[mi];
	}

	double V1Len = 0;
	double V2Len = 0;
	for (int mi = 0; mi < m; ++mi) {
		V1Len += V1[mi] * V1[mi];
		V2Len += V2[mi] * V2[mi];
	}
	V1Len = sqrt(V1Len);
	V2Len = sqrt(V2Len);

	if (V1Len == 0.0 || V2Len == 0.0) return ILLEGAL_COS;
	else return V1V2Multiply / V1Len / V2Len;
}


double calculateCos_PointToVector(const vector<double> &point, const vector<double> &refPoint, const vector<double> &vec) {

	int m = vec.size();
	if (m != point.size() || m != refPoint.size()) return ILLEGAL_COS;

	vector<double> V1(m);
	double V1V2Multiply = 0;
	for (int mi = 0; mi < m; ++mi) {
		V1[mi] = point[mi] - refPoint[mi];
		V1V2Multiply += V1[mi] * vec[mi];
	}

	double V1Len = 0;
	double V2Len = 0;
	for (int mi = 0; mi < m; ++mi) {
		V1Len += V1[mi] * V1[mi];
		V2Len += vec[mi] * vec[mi];
	}
	V1Len = sqrt(V1Len);
	V2Len = sqrt(V2Len);

	if (V1Len == 0.0 || V2Len == 0.0) return ILLEGAL_COS;
	else return V1V2Multiply / V1Len / V2Len;
}


bool judgeInBound_M2M(solution &sol, const vector<double> &refPoint, const vector<double> &myWV, const vector< vector<double> > &otherWV) {
	if (sol.fitness.size() != refPoint.size()) return false;
	
	sol.cosAngleToWV = calculateCos_PointToVector(sol.fitness, refPoint, myWV);
	if (ILLEGAL_COS == sol.cosAngleToWV) {
		sol.inBoundFlag = false;
		return sol.inBoundFlag;
	}

	sol.inBoundFlag = true;
	for (int otherIndex = 0; otherIndex < otherWV.size(); ++otherIndex) {
		double otherCos = calculateCos_PointToVector(sol.fitness, refPoint, otherWV[otherIndex]);
		if (ILLEGAL_COS == otherCos) {
			sol.inBoundFlag = false;
			return sol.inBoundFlag;
		}

		if (otherCos > sol.cosAngleToWV) {
			sol.inBoundFlag = false;
			break;
		}
	}

	return sol.inBoundFlag;
}


double scalarFunc_Chebyshev(const vector<double> &fitness, const vector<double> &refPoint, const vector<double> &weiVector)
{
	int m = weiVector.size();
	double minWeiFitness = 1E300;
	double theWeiFitness;
	for (int mi = 0; mi < m; ++mi) {
		if (weiVector[mi] != 0) {
			theWeiFitness = 1.0 / weiVector[mi] * (fitness[mi] - refPoint[mi]);
			if (theWeiFitness < minWeiFitness) {
				minWeiFitness = theWeiFitness;
			}
		}
	}
	return minWeiFitness;
}


bool judgeFitCorrect(const problem &inst, solution &sol) {
	int m = sol.fitness.size();
	vector<double> oriFit = sol.fitness;
	sol.calculateFitness(inst);

	bool correctFlag = true;
	for (int mi = 0; mi < m; mi++) {
		if (sol.fitness[mi] != oriFit[mi]) {
			correctFlag = false;
			break;
		}
	}
	return correctFlag;
}

bool acceptanceCriterion1(const solution &sol, const archive &archiveA) {
	//Acceptance criterion 1: if a solution has a better scalar function value than all solutions in archive, then accept it. 
	//Out-bound is not alowed, except when no solution in the archive is in-bound.
	bool acceptFlag = false;
	if (sol.inBoundFlag || 0 == archiveA.inBoundSolNum) {
		if (sol.scalarFunctionValue > archiveA.bestScalarFuncValue) {
			acceptFlag = true;
		}
	}
	return acceptFlag;
}

bool acceptanceCriterion2(const solution &sol, const archive &archiveA) {
	//Acceptance criterion 2: if a solution is not dominated by any solutions in archive, then accept it. 
	//Repeatition is not alowed.
	//Out-bound is not alowed, except when no solution in the archive is in-bound.
	bool acceptFlag = false;
	if (true == sol.inBoundFlag || 0 == archiveA.inBoundSolNum) {
		if (!archiveA.judgeBeDomd_MaxCase(sol) && !archiveA.judgeRepeatition(sol)) {
			acceptFlag = true;
		}
	}
	return acceptFlag;
}

bool acceptanceCriterion1_AlowOutBound(const solution &sol, const archive &archiveA) {
	bool acceptFlag = false;
	//If a solution has a better scalar function value than all solutions in archive, 
	//then no matter whether it is in-bound or out-bound, accept it. 
	if (sol.scalarFunctionValue > archiveA.bestScalarFuncValue) {
		acceptFlag = true;
	}

	return acceptFlag;
}

solution genNeighborSol(const problem &inst, const solution &sol, const vector< vector <int> > &twoOptFlipCityList, const int k) {
	int P = twoOptFlipCityList[k][0];
	int Q = twoOptFlipCityList[k][1];

	solution solPrime = sol;
	for (int mi = 0; mi < inst.m; mi++) {
		solPrime.fitness[mi] += moveCost2Opt(mi, P, Q, sol, inst);
	}
	move2Opt(P, Q, solPrime);
	return solPrime;
}

solution genNeighborSol_WithoutFitness(const problem &inst, const solution &sol, const vector< vector <int> > &twoOptFlipCityList, const int k) {
	int P = twoOptFlipCityList[k][0];
	int Q = twoOptFlipCityList[k][1];

	solution solPrime = sol;
	move2Opt(P, Q, solPrime);
	return solPrime;
}

bool getMyProcsInfo(const char* wvFilename, const int numProcs, const int myId, const vector<double> &refPoint, My_Procs_Info &info) {
	//load weight vectors from file abd store in myProcsInfo
	int m = refPoint.size();
	info.numProcs = numProcs;
	info.weightVectorFilename = wvFilename;
	info.myId = myId;
	info.myWeightVector.resize(m); //my own weight vector
	info.otherId.resize(numProcs - 1); //all other processes' IDs, sorted based on the distance (angle) to my own weight vector
	info.otherWeightVector.resize(numProcs - 1); //all other processes' weight vectors, sorted based on the distance (angle) to my own weight vector
	for (int index = 0; index < numProcs - 1; index++) {
		info.otherWeightVector[index].resize(m);
	}
	info.refPoint = refPoint;

	FILE *wvFile = fopen(wvFilename, "r");
	if (wvFile == NULL) {
		printf("ERROR: cannot open weight vector file %s\n", wvFilename);
		return false;
	}

	int IdReadBuffer;
	int otherIndex = 0;
	for (int line = 0; line < numProcs; line++) {
		if (fscanf(wvFile, "%d", &IdReadBuffer) != 1) {
			printf("ERROR: wrong format in weight vector file %s\n", wvFilename);
			fclose(wvFile);
			return false;
		}
		if (IdReadBuffer == myId) {
			for (int mi = 0; mi < m; ++mi) {
				if (fscanf(wvFile, "%lf", &(info.myWeightVector[mi])) != 1) {
					printf("ERROR: wrong format in weight vector file %s\n", wvFilename);
					fclose(wvFile);
					return false;
				}
			}
		}
		else {
			info.otherId[otherIndex] = IdReadBuffer;
			for (int mi = 0; mi < m; ++mi) {
				if (fscanf(wvFile, "%lf", &(info.otherWeightVector[otherIndex][mi])) != 1) {
					printf("ERROR: wrong format in weight vector file %s\n", wvFilename);
					fclose(wvFile);
					return 1;
				}
			}
			otherIndex++;
		}
	}
	fclose(wvFile);

	//sort 'otherId' and 'otherWeightVector' based on its distance (angle) to 'myWeightVector'
	multimap<double, vector<double> > sortSpace;
	for (int index = 0; index < numProcs - 1; index++) {
		double tempCos = calculateCos_VectorToVector(info.otherWeightVector[index], info.myWeightVector);

		vector<double> temp(m + 1);
		temp[0] = info.otherId[index];
		for (int mi = 0; mi < m; mi++) {
			temp[1 + mi] = info.otherWeightVector[index][mi];
		}

		sortSpace.insert(pair< double, vector<double> >(tempCos, temp));
	}
	otherIndex = 0;
	for (multimap<double, vector<double> >::reverse_iterator it = sortSpace.rbegin(); it != sortSpace.rend(); it++) {
		info.otherId[otherIndex] = it->second[0];
		for (int mi = 0; mi < m; mi++) {
			info.otherWeightVector[otherIndex][mi] = it->second[1 + mi];
		}
		otherIndex++;
	}

	return true;
}


