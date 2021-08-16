#include "pplsd.h"

archive PPLSD(const problem &inst, const char* wvFilename, const archive &initialArchive, const int maxRuntime, const int numProcs, const int myId, const vector<double> refPoint) {

	int m = inst.m;
	int n = inst.n;

	//generate the 2-Opt flipping city list
	int neighborCount = 0;
	vector< vector<int> > twoOptFlipCityList;
	for (int P = 0; P < n; P++) {
		for (int Q = P + 2; Q < n; Q++) {
			vector<int> tempPQ(2);
			tempPQ[0] = P;
			tempPQ[1] = Q;
			twoOptFlipCityList.push_back(tempPQ);
			neighborCount++;
		}
	}
	int neighborhoodSize = neighborCount;
	vector< vector<double> > neighborhoodRecord(neighborhoodSize);//the room to store neighboring solutions' fitness
	archive emptyArchive;

	// read weight vector file and get all information about my process
	My_Procs_Info myInfo;
	if (!getMyProcsInfo(wvFilename, numProcs, myId, refPoint, myInfo)) {
		return emptyArchive;
	}

	//initialization archive
	archive archiveA;
	for (list<solution>::const_iterator itSol = initialArchive.solList.begin(); itSol != initialArchive.solList.end(); itSol++) {
		solution iniSol = *itSol;
		iniSol.scalarFunctionValue = scalarFunc_Chebyshev(iniSol.fitness, refPoint, myInfo.myWeightVector);
		iniSol.inBoundFlag = judgeInBound_M2M(iniSol, refPoint, myInfo.myWeightVector, myInfo.otherWeightVector);
		iniSol.explored = false;
		archiveA.updateArchive(iniSol);
	}

	double functionEvaluateCount = 0;
	double currentTime = 0.0;
	bool timeUpFlag = false;
	start_timer();

	while (!timeUpFlag) {

		archiveA.countUnexploredSolNum();
		if (!archiveA.solList.empty() && 0 == archiveA.unexploredSolNum) {
			break;
		}

		//select a unexplored solution from archive with the max scalar function value as the current solution 'sol'
		double unexploredBestValue = -DBL_MAX;
		list<solution>::iterator itSelectedSol;
		for (list<solution>::iterator itSol = archiveA.solList.begin(); itSol != archiveA.solList.end(); itSol++) {
			if (!itSol->explored) {
				if (itSol->scalarFunctionValue > unexploredBestValue) {
					unexploredBestValue = itSol->scalarFunctionValue;
					itSelectedSol = itSol;
				}
			}
		}
		solution sol = *itSelectedSol;
		bool surviveFlag = true;

		//first round neighborhood exploration using acceptance criterion 1 and first-improve strategy
		//traverse the neighboring solution of current solution 'sol'
		bool successFlag = false;
		for (int k = 0; k < neighborhoodSize; k++) {

#if defined(__linux__) || defined(__linux)
			currentTime = return_CPU_time();
#else
			currentTime = return_elapsed_time();
#endif
			if (currentTime >= maxRuntime) {
				timeUpFlag = true;
				break;
			}

			//generate the k-th neigboring solution 'solPrime'
			solution solPrime = genNeighborSol(inst, sol, twoOptFlipCityList, k);
			solPrime.scalarFunctionValue = scalarFunc_Chebyshev(solPrime.fitness, refPoint, myInfo.myWeightVector);
			solPrime.inBoundFlag = judgeInBound_M2M(solPrime, refPoint, myInfo.myWeightVector, myInfo.otherWeightVector);

			functionEvaluateCount++;

			//storage the fitness of the k-th neighboring solution
			neighborhoodRecord[k] = solPrime.fitness;

			//If the neighboring solution 'solPrime' is dominated by the original solution 'sol', skip it
			if (sol.judgeDominate(solPrime)) {
				continue;
			}

			//judge whether 'solPrime' satisfies the acceptance criterion 1
			if (acceptanceCriterion1(solPrime, archiveA)) {

				//update archive by 'solPrime'
				solPrime.explored = false;
				archiveA.updateArchive(solPrime);

				//if the current solution 'sol' is dominated by the newly added 'solPrime', 
				//then adding 'solPrine' to archive will remove 'sol', in other words, 'sol' does not survive.
				if (solPrime.judgeDominate(sol)) {
					surviveFlag = false;
				}

				//Acceptance criterion 1 applies first-improve strategy
				//the first neighboring solution that statisfy acceptance criterion 1 has been found , stop avluating the remain neighboring solutuons
				successFlag = true;
				break;

			}
		}


		//If first round fails, start second rpund neighborhood exploration using acceptance criterion 2 and best-improve strategy
		if (!successFlag) {

			//traverse the neighboring solution of current solution 'sol'
			for (int k = 0; k < neighborhoodSize; k++) {

#if defined(__linux__) || defined(__linux)
				currentTime = return_CPU_time();
#else
				currentTime = return_elapsed_time();
#endif
				if (currentTime >= maxRuntime) {
					timeUpFlag = true;
					break;
				}

				//generate the k-th neigboring solution 'solPrime'
				//since the fitness of the k-th neighboring solution is already in the storage, we generate 'solPrime' without function evaluation
				solution solPrime = genNeighborSol_WithoutFitness(inst, sol, twoOptFlipCityList, k);
				solPrime.fitness = neighborhoodRecord[k];
				solPrime.scalarFunctionValue = scalarFunc_Chebyshev(solPrime.fitness, refPoint, myInfo.myWeightVector);
				solPrime.inBoundFlag = judgeInBound_M2M(solPrime, refPoint, myInfo.myWeightVector, myInfo.otherWeightVector);

				//If the neighboring solution 'solPrime' is dominated by the original solution 'sol', skip it
				if (sol.judgeDominate(solPrime)) {
					continue;
				}

				//judge whether 'solPrime' satisfies the acceptance criterion 1
				if (acceptanceCriterion2(solPrime, archiveA)) {

					//update archive by 'solPrime'
					solPrime.explored = false;
					archiveA.updateArchive(solPrime);

					//if the current solution 'sol' is dominated by the newly added 'solPrime', 
					//then adding 'solPrine' to archive will remove 'sol', in other words, 'sol' does not survive.
					if (solPrime.judgeDominate(sol)) {
						surviveFlag = false;
					}

					//Do NOT break the neighborhood traverse for loop,
					//since acceptance criterion 2 applies best-improve strategy,
					//in other words, all the neighboring solution that statisfy acceptance criterion 2 will be accepted to archive
				}
			}
		}

		//After the neighborhood exploration, if the current solution 'sol' survives, which means it is still in archive, then mark it as explored
		if (surviveFlag) {
			itSelectedSol->explored = true;
		}
	}

	return archiveA;
}









