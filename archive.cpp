#include "archive.h"

using namespace std;

archive::archive() {
	unexploredSolNum = 0;
	inBoundSolNum = 0;
	bestScalarFuncValue = -DBL_MAX;
}

void archive::addToArchive(const solution &theSol) {
	solList.push_back(theSol);
	return;
}

void archive::deleteFromArchive(list<solution>::iterator &solIt) {
	solList.erase(solIt);
	return;
}

bool archive::judgeRepeatition(const solution &theSol) const {
	if (solList.empty()) {
		return false;
	}

	for (list<solution>::const_iterator it = solList.begin(); it != solList.end(); it++) {
		if (it->judgeSameSol(theSol)) {
			return true;
		}
	}
	return false;
}

bool archive::judgeBeDomd_MaxCase(const solution &theSol) const {
	if (solList.empty()) {
		return false;
	}

	for (list<solution>::const_iterator it = solList.begin(); it != solList.end(); it++) {
		if (it->judgeDominate(theSol)) {
			return true;
		}
	}
	return false;
}

void archive::updateArchive(const solution &theSol) {

	//delete the solutions that are dominated by theSol

	list<list<solution>::iterator> toBeDelete;
	bool bestChangeFlag = false;
	for (list<solution>::iterator it = solList.begin(); it != solList.end(); it++) {
		if (theSol.judgeDominate(*it)) {
			toBeDelete.push_back(it);
		}
	}
	for (list<list<solution>::iterator>::iterator it2 = toBeDelete.begin(); it2 != toBeDelete.end(); it2++) {
		if ((*it2)->inBoundFlag) {
			inBoundSolNum--;
		}
		if ((*it2)->scalarFunctionValue == bestScalarFuncValue) {
			bestChangeFlag = true;
		}
		deleteFromArchive(*it2);
	}
	if (bestChangeFlag) {
		getBestScalarFuncValue();
	}
	
	// add solution to archive
	addToArchive(theSol);
	if (theSol.inBoundFlag) {
		inBoundSolNum++;
	}
	if (theSol.scalarFunctionValue > bestScalarFuncValue) {
		bestScalarFuncValue = theSol.scalarFunctionValue;
	}
	return;
}

int archive::countUnexploredSolNum() {
	unexploredSolNum = 0;
	for (list<solution>::iterator it = solList.begin(); it != solList.end(); it++) {
		if (!(it->explored)) {
			unexploredSolNum++;
		}
	}
	return unexploredSolNum;
}

int archive::countInBoundSolNum() {
	inBoundSolNum = 0;
	for (list<solution>::iterator it = solList.begin(); it != solList.end(); it++) {
		if (it->inBoundFlag) {
			inBoundSolNum++;
		}
	}
	return inBoundSolNum;
}

double archive::getBestScalarFuncValue() {
	bestScalarFuncValue = -DBL_MAX;
	for (list<solution>::iterator it = solList.begin(); it != solList.end(); it++) {
		if (it->scalarFunctionValue > bestScalarFuncValue) {
			bestScalarFuncValue = it->scalarFunctionValue;
		}
	}
	return bestScalarFuncValue;
}

bool archive::loadFromFile_TransferToMaxProblem(const char* filename) {
	//load the archive from file and transfer the archive to a maximization problem's archive by multiplying -1

	FILE *iniAFile = fopen(filename, "r");
	if (iniAFile == NULL) {
		printf("ERROR: cannot open file %s\n", filename);
		return false;
	}

	solList.clear();
	
	int archiveSize = 0;
	if (fscanf(iniAFile, "%d", &archiveSize) == 0) {
		printf("ERROR: cannot read the archive size in %s\n", filename);
		fclose(iniAFile);
		return false;
	}

	if (archiveSize > 0) {
		int m;
		if (fscanf(iniAFile, "%d", &m) == 0) {
			printf("ERROR: cannot read the objective number m in %s\n", filename);
			fclose(iniAFile);
			return false;
		}

		int n;
		if (fscanf(iniAFile, "%d", &n) == 0) {
			printf("ERROR: cannot read the solution dimension n in %s\n", filename);
			fclose(iniAFile);
			return false;
		}


		solution iniSol;
		double tempFitness;
		iniSol.fitness.resize(m);
		iniSol.x.resize(n);

		for (int line = 0; line < archiveSize; line++) {
			for (int mi = 0; mi < m; mi++) {
				if (fscanf(iniAFile, "%lf", &tempFitness) == 0) {
					printf("ERROR: cannot read the %d-th solution in %s\n", line, filename);
					fclose(iniAFile);
					return false;
				}
				iniSol.fitness[mi] = -1 * tempFitness;//transfer to a maximization problem by multiply -1
			}

			for (int i = 0; i < n; i++) {
				if (fscanf(iniAFile, "%d", &iniSol.x[i]) == 0) {
					printf("ERROR: cannot read the %d-th solution in %s\n", line, filename);
					fclose(iniAFile);
					return false;
				}
			}

			updateArchive(iniSol);
		}
	}
	
	fclose(iniAFile);
	return true;
}

bool archive::printToFile_RestoreMinProblem(const char* filename) {
	//print the archive to file and restore the original minimization problem

	FILE *resultSolFile = fopen(filename, "w");
	if (resultSolFile == NULL) {
		printf("ERROR: cannot open file %s\n", filename);
		return false;
	}

	fprintf(resultSolFile, "%d\n", solList.size());
	if (!solList.empty()) {
		int m = solList.begin()->fitness.size();
		int n = solList.begin()->x.size();
		fprintf(resultSolFile, "%d %d\n", m, n);
		for (list<solution>::iterator itSol = solList.begin(); itSol != solList.end(); itSol++) {
			for (int mi = 0; mi < m; mi++) {
				fprintf(resultSolFile, "%.0f ", -1 * itSol->fitness[mi]);//restore the original minimization problem
			}
			for (int i = 0; i < n; ++i) {
				fprintf(resultSolFile, " %d", itSol->x[i]);
			}
			fprintf(resultSolFile, "\n");
		}
	}
	fclose(resultSolFile);
	return true;
}


