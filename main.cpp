
#include "problem.h"
#include "solution.h"
#include "archive.h"
#include "misc.h"
#include "pplsd.h"

#include "mpi.h"

using namespace std;

int main(int argc, char *argv[])
{
	int myId, numProcs, hostNameLen;
	char myHostName[MPI_MAX_PROCESSOR_NAME];
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myId);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Get_processor_name(myHostName, &hostNameLen);

	if (argc < 4) {
		printf("ERROR: Arguments missing! The arguments should be <problem_filename> <weight_vector_filename> <max_runtime>\n");
		MPI_Finalize();
		return 1;
	}
	//argv[1] is the problem instance filename
	//argv[2] is the weight vector filename
	int maxRuntime = strtol(argv[3], NULL, 10);

	problem inst;
	if (!inst.readProblem_TransferToMaxProblem(argv[1])) {
		printf("ERROR: cannot open problem instance file %s\n", argv[1]);
		MPI_Finalize();
		return 1;
	}

	int n = inst.n;
	int m = inst.m;

	//the initial archive contains only one randomly generated solution
	srand(1);//use fixed random seed to make sure that different processes have the same initial solution, hence the same reference point  
	solution iniSol;
	iniSol.randomInitialize(inst);
	archive iniArchive;
	iniArchive.updateArchive(iniSol);

	//the reference point is set to be equal to the function vector of the randomly generated solution
	vector<double> refPoint(m);
	for (int mi = 0; mi < m; mi++) {
		refPoint[mi] = iniSol.fitness[mi];
	}

	archive result = PPLSD(inst, argv[2], iniArchive, maxRuntime, numProcs, myId, refPoint);

	char resultFilename[200];
	sprintf(resultFilename, "result_process%d.txt", myId);
	result.printToFile_RestoreMinProblem(resultFilename);

	printf("Process %d finished, archive size = %d, output file : result_process%d.txt\n", myId, result.solList.size(), myId);

	MPI_Barrier(MPI_COMM_WORLD);

	//the process 0 will combine the results of all processes
	if (myId == 0) {
		for (int procsId = 1; procsId < numProcs; procsId++) {
			char otherFilename[200];
			sprintf(otherFilename, "result_process%d.txt", procsId);

			archive tempArchive;
			tempArchive.loadFromFile_TransferToMaxProblem(otherFilename);

			for (list<solution>::iterator itSol = tempArchive.solList.begin(); itSol != tempArchive.solList.end(); itSol++) {
				if (!result.judgeBeDomd_MaxCase(*itSol)) {
					result.updateArchive(*itSol);
				}
			}
		}

		char finalResultFilename[200];
		sprintf(finalResultFilename, "result_final.txt");
		result.printToFile_RestoreMinProblem(finalResultFilename);

		printf("All processes finished, final archive size = %d, output file : result_final.txt\n", result.solList.size());
	}

	MPI_Finalize();

	return 0;
}








