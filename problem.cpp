#include "problem.h"

using namespace std;

double dtrunc(double x)
{
	int k;

	k = (int)x;
	x = (double)k;
	return x;
}

bool problem::readProblem_TransferToMaxProblem(const char* filename)
{
	ifstream tspStream(((string)filename).c_str());
	if (!tspStream) {
		return false;
	}
	char key[256];
	int noUse;

	while (strncmp(key, "TYPE", 4)) {
		tspStream >> key;
	}
	memset(key, 0, sizeof(key));
	tspStream >> key;
	if (!strncmp(key, ":", 1)) {
		tspStream >> key;
	}
	if (strncmp(key, "mTSP", 4)) {
		return false;
	}

	while (strncmp(key, "OBJECTIVE_NUM", 13)) {
		tspStream >> key;
	}
	memset(key, 0, sizeof(key));
	tspStream >> key;
	if (!strncmp(key, ":", 1)) {
		tspStream >> key;
	}
	m = atoi(key);

	while (strncmp(key, "DIMENSION", 9)) {
		tspStream >> key;
	}
	memset(key, 0, sizeof(key));
	tspStream >> key;
	if (!strncmp(key, ":", 1)) {
		tspStream >> key;
	}
	n = atoi(key);

	while (strncmp(key, "EDGE_WEIGHT_TYPE", 16)) {
		tspStream >> key;
	}
	memset(key, 0, sizeof(key));
	tspStream >> key;
	if (!strncmp(key, ":", 1)) {
		tspStream >> key;
	}

	if (!strncmp(key, "EUC_2D", 6)) {
		edgeWeightType = CC_EUCLIDEAN;
	}
	else if (!strncmp(key, "CEIL_2D", 7)) {
		edgeWeightType = CC_EUCLIDEAN_CEIL;
	}
	else if (!strncmp(key, "GEO", 3)) {
		edgeWeightType = CC_GEOGRAPHIC;
	}
	else if (!strncmp(key, "ATT", 3)) {
		edgeWeightType = CC_ATT;
	}
	else {
		return false;
	}

	dist.resize(m);
	for (int mi = 0; mi < m; ++mi) {
		dist[mi].resize(n);
		for (int i = 0; i < n; ++i) {
			dist[mi][i].resize(n);
		}
	}

	//if (Norm != CC_MATRIXNORM)
	while (strncmp(key, "NODE_COORD_SECTION", 18)) {
		tspStream >> key;
	}
	X.resize(m);
	Y.resize(m);
	for (int mi = 0; mi < m; ++mi) {
		X[mi].resize(n);
		Y[mi].resize(n);
	}

	for (int i = 0; i < n; i++) {
		tspStream >> noUse;
		for (int mi = 0; mi < m; mi++) {
			tspStream >> X[mi][i] >> Y[mi][i];
		}
	}

	switch (edgeWeightType)
	{
	case CC_EUCLIDEAN:
	{
		for (int mi = 0; mi < m; mi++) {
			for (int i = 0; i < n; i++) {
				dist[mi][i][i] = 0;
				for (int j = 0; j < i; j++) {
					dist[mi][i][j] = (int)floor(sqrt((X[mi][i] - X[mi][j])*(X[mi][i] - X[mi][j]) + (Y[mi][i] - Y[mi][j])*(Y[mi][i] - Y[mi][j])) + 0.5);
					dist[mi][j][i] = dist[mi][i][j];
				}
			}
		}
		break;
	}
	case CC_EUCLIDEAN_CEIL:
	{
		for (int mi = 0; mi < m; mi++) {
			for (int i = 0; i < n; i++) {
				dist[mi][i][i] = 0;
				for (int j = 0; j < i; j++) {
					dist[mi][i][j] = (int)ceil(sqrt((X[mi][i] - X[mi][j])*(X[mi][i] - X[mi][j]) + (Y[mi][i] - Y[mi][j])*(Y[mi][i] - Y[mi][j])));
					dist[mi][j][i] = dist[mi][i][j];
				}
			}
		}
		break;
	}
	case CC_GEOGRAPHIC:
	{
		for (int mi = 0; mi < m; mi++) {
			for (int i = 0; i < n; i++) {
				dist[mi][i][i] = 0;
				for (int j = 0; j < i; j++) {
					double deg, min;
					double lati, latj, longi, longj;
					double q1, q2, q3;
					double x1 = X[mi][i], x2 = X[mi][j], yy1 = Y[mi][i], yy2 = Y[mi][j];

					deg = dtrunc(x1);
					min = x1 - deg;
					lati = GH_PI * (deg + 5.0 * min / 3.0) / 180.0;
					deg = dtrunc(x2);
					min = x2 - deg;
					latj = GH_PI * (deg + 5.0 * min / 3.0) / 180.0;

					deg = dtrunc(yy1);
					min = yy1 - deg;
					longi = GH_PI * (deg + 5.0 * min / 3.0) / 180.0;
					deg = dtrunc(yy2);
					min = yy2 - deg;
					longj = GH_PI * (deg + 5.0 * min / 3.0) / 180.0;

					q1 = cos(longi - longj);
					q2 = cos(lati - latj);
					q3 = cos(lati + latj);
					dist[mi][i][j] = (int)(6378.388 * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
					dist[mi][j][i] = dist[mi][i][j];
				}
			}
		}
		break;
	}
	case CC_ATT:
	{
		for (int mi = 0; mi < m; mi++) {
			for (int i = 0; i < n; i++) {
				dist[mi][i][i] = 0;
				for (int j = 0; j < i; j++) {
					double xd = X[mi][i] - X[mi][j];
					double yd = Y[mi][i] - Y[mi][j];
					double rij = sqrt((xd * xd + yd * yd) / 10.0);
					double tij = dtrunc(rij);

					if (tij < rij)
						dist[mi][i][j] = (int)tij + 1;
					else
						dist[mi][i][j] = (int)tij;
					dist[mi][j][i] = dist[mi][i][j];
				}
			}
		}
		break;
	}
	}

	//since TSP is a minimization problem, we transfer it to a maximization problem by multipling -1
	for (int mi = 0; mi < m; mi++) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				dist[mi][i][j] = -1 * dist[mi][i][j];
			}
		}
	}

	return true;
}


