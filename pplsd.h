#pragma once

#include "problem.h"
#include "solution.h"
#include "archive.h"
#include "misc.h"

using namespace std;

archive PPLSD(const problem &inst, const char* wvFilename, const archive &initialArchive, const int maxRuntime, const int numProcs, const int myId, const vector<double> refPoint);






