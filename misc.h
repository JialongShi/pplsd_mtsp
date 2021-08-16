#pragma once

#include "problem.h"
#include "solution.h"
#include "archive.h"

#include <map>
#include <limits.h>
#include <ctime>
#if defined(__linux__) || defined(__linux)
#include <sys/time.h>
#include <sys/times.h>
#include <unistd.h>
#else
#include <windows.h>
#endif

using namespace std;

#if defined(__linux__) || defined(__linux)
static struct tms  Old_Time;
static struct tms  New_Time;
static int ClkTck;
static clock_t clock_t_start;
static clock_t clock_t_end;

void start_timer();

double return_CPU_time();

double return_elapsed_time();

#else
static time_t Old_Time;
static time_t New_Time;

void start_timer();

double return_elapsed_time();

#endif




double moveCost2Opt(const int &mi, const int &P, const int &Q, const solution &sol, const problem &inst);

bool move2Opt(const int &P, const int &Q, solution &sol);

double calculateCos_VectorToVector(const vector<double> &V1, const vector<double> &V2);

double calculateCos_PointToVector(const vector<double> &point, const vector<double> &refPoint, const vector<double> &vec);

bool judgeInBound_M2M(solution &sol, const vector<double> &refPoint, const vector<double> &myWV, const vector< vector<double> > &otherWV);

double scalarFunc_Chebyshev(const vector<double> &fitness, const vector<double> &refPoint, const vector<double> &weiVector);

bool judgeFitCorrect(const problem &inst, solution &sol);

bool acceptanceCriterion1(const solution &sol, const archive &archiveA);

bool acceptanceCriterion2(const solution &sol, const archive &archiveA);

bool acceptanceCriterion1_AlowOutBound(const solution &sol, const archive &archiveA);

solution genNeighborSol(const problem &inst, const solution &sol, const vector< vector <int> > &twoOptFlipCityList, const int k);

solution genNeighborSol_WithoutFitness(const problem &inst, const solution &sol, const vector< vector <int> > &twoOptFlipCityList, const int k);



struct My_Procs_Info {
	int numProcs;
	string weightVectorFilename;
	int myId;
	vector<double> myWeightVector;
	vector<int> otherId;
	vector< vector<double> > otherWeightVector;
	vector<double> refPoint;
};

bool getMyProcsInfo(const char* wvFilename, const int numProcs, const int myId, const vector<double> &refPoint, My_Procs_Info &info);


