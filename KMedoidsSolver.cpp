#include "KMedoidsSolver.h"

using namespace std;

const double MAX_FLOAT = std::numeric_limits<double>::max();

KMedoidsSolver::KMedoidsSolver(DataFrame dataFrame, int* solution)
: KCenterSolver(dataFrame, solution) {
	createCenters();
	calculateCost();
}

KMedoidsSolver::~KMedoidsSolver() {

}

std::vector<int>* KMedoidsSolver::getClusters() const {
	return this->clusters;
}

double KMedoidsSolver::getMedian(std::vector<int> cl, int j) {
	const int clusterSize = cl.size();
	double values[clusterSize];
	double** data = this->dataFrame.getData();
	
	for(int i = 0; i < clusterSize; i++) {
		values[i] = data[cl[i]][j];
	}

	int pos = clusterSize/2;

	nth_element(values, values + pos, values + clusterSize);
	double median = values[pos];
	
	return median;
}

void KMedoidsSolver::updateClusters(int p, int c1, int c2) {
	int y = 0;
	
	for(int i = 0; i < this->clusters[c1].size(); i++) {
		if(this->clusters[c1][i] == p) {
			y = i;
		}
	}

	this->clusters[c2].push_back(p);
	this->clusters[c1].erase(this->clusters[c1].begin()+y);
}

void KMedoidsSolver::updateClusters(int p1, int p2, int c1, int c2) {
	int y1 = 0;
	int y2 = 0;
	
	for(int i = 0; i < this->clusters[c1].size(); i++) {
		if(this->clusters[c1][i] == p1) {
			y1 = i;
		}
	}

	for(int i = 0; i < this->clusters[c2].size(); i++) {
		if(this->clusters[c2][i] == p2) {
			y2 = i;
		}
	}

	this->clusters[c1][y1] = p2;
	this->clusters[c2][y2] = p1;
}

void KMedoidsSolver::createCenters() {
	int N = this->dataFrame.getInstance().N;
	int M = this->dataFrame.getInstance().M;
	int D = this->dataFrame.getInstance().D;
	double** data = this->dataFrame.getData();

	int centers[M];
	this->centroid = new double*[M];
	this->clusters = new vector<int>[M];
	this->cardinality = new int[M];

	double min;
	double dist;

	for(int i = 0; i < M; i++) {
		this->centroid[i] = new double[D];
		this->cardinality[i] = 0;
	}

	for(int i = 0; i < N; i++) {
		this->cardinality[solution[i]]++;
		this->clusters[solution[i]].push_back(i);
	}

	for(int i = 0; i < M; i++) {
		for(int j = 0; j < D; j++) {
			this->centroid[i][j] = getMedian(this->clusters[i], j);
		}
	}
	
	for(int i = 0; i < M; i++) {
		min = MAX_FLOAT;
		// setting the real centroid as the closest point to the artificial centroid calculated above
		for(int k = 0; k < this->clusters[i].size(); k++) {
			dist = getDistance(data[this->clusters[i][k]], centroid[i], D);
			if(dist < min) {
				min = dist;
				centers[i] = this->clusters[i][k];
			}
		}
	}

	for(int i = 0; i < M; i++) {
		for(int j = 0; j < D; j++) {
			this->centroid[i][j] = data[centers[i]][j];
		}
	}
}

void KMedoidsSolver::calculateNewCentroids(int p, int c2, double* newCentroid1, double* newCentroid2) {
	int d = this->dataFrame.getInstance().D;
	int n = this->dataFrame.getInstance().N;
	int c1 = this->solution[p];

	double** data = this->dataFrame.getData();
	int y = 0;

	for(int i = 0; i < this->clusters[c1].size(); i++) {
		if(this->clusters[c1][i] == p) {
			y = i;
		}
	}

	vector<int> cl1 = this->clusters[c1];
	vector<int> cl2 = this->clusters[c2];

	cl2.push_back(p);
	cl1.erase(cl1.begin()+y);

	for(int j = 0; j < d; j++) {
		newCentroid1[j] = getMedian(cl1, j);
		newCentroid2[j] = getMedian(cl2, j);
	}

	double dist;
	int center1;
	int center2;

	double min = MAX_FLOAT;
	for(int k = 0; k < cl1.size(); k++) {
		dist = getDistance(data[cl1[k]], newCentroid1, d);
		if(dist < min) {
			min = dist;
			center1 = cl1[k];
		}
	}

	min = MAX_FLOAT;
	for(int k = 0; k < cl2.size(); k++) {
		dist = getDistance(data[cl2[k]], newCentroid2, d);
		if(dist < min) {
			min = dist;
			center2 = cl2[k];
		}
	}

	for(int j = 0; j < d; j++) {
		newCentroid1[j] = data[center1][j];
		newCentroid2[j] = data[center2][j];
	}
}

void KMedoidsSolver::relocate(int p, int c2) {
	int N = this->dataFrame.getInstance().N;
	int M = this->dataFrame.getInstance().M;
	int D = this->dataFrame.getInstance().D;
	int c1 = this->solution[p];
	double** data = this->dataFrame.getData();

	this->cardinality[c1] = this->cardinality[c1] - 1;
	this->cardinality[c2] = this->cardinality[c2] + 1;

	updateClusters(p, c1, c2);

	for(int j = 0; j < D; j++) {
		this->centroid[c1][j] = getMedian(this->clusters[c1], j);
		this->centroid[c2][j] = getMedian(this->clusters[c2], j);
	}

	double min;
	double dist;

	int centers[M];

	for(int i = 0; i < M; i++) {
		min = MAX_FLOAT;
		// setting the real centroid as the closest point to the artificial centroid calculated above
		for(int k = 0; k < this->clusters[i].size(); k++) {
			dist = getDistance(data[this->clusters[i][k]], this->centroid[i], D);
			if(dist < min) {
				min = dist;
				centers[i] = this->clusters[i][k];
			}
		}
	}

	for(int i = 0; i < M; i++) {
		for(int j = 0; j < D; j++) {
			this->centroid[i][j] = data[centers[i]][j];
		}
	}
}

void KMedoidsSolver::calculateNewCentroidsSwap(int p1, int p2, double* newCentroid1, double* newCentroid2) {
	int d = this->dataFrame.getInstance().D;
	int n = this->dataFrame.getInstance().N;
	int c1 = this->solution[p1];
	int c2 = this->solution[p2];

	double** data = this->dataFrame.getData();
	int y1 = 0;
	int y2 = 0;

	for(int i = 0; i < this->clusters[c1].size(); i++) {
		if(this->clusters[c1][i] == p1) {
			y1 = i;
		}
	}

	for(int i = 0; i < this->clusters[c2].size(); i++) {
		if(this->clusters[c2][i] == p2) {
			y2 = i;
		}
	}

	vector<int> cl1 = this->clusters[c1];
	vector<int> cl2 = this->clusters[c2];

	cl1[y1] = p2;
	cl2[y2] = p1;

	for(int j = 0; j < d; j++) {
		newCentroid1[j] = getMedian(cl1, j);
		newCentroid2[j] = getMedian(cl2, j);
	}

	double dist;
	int center1;
	int center2;

	double min = MAX_FLOAT;
	for(int k = 0; k < cl1.size(); k++) {
		dist = getDistance(data[cl1[k]], newCentroid1, d);
		if(dist < min) {
			min = dist;
			center1 = cl1[k];
		}
	}

	min = MAX_FLOAT;
	for(int k = 0; k < cl2.size(); k++) {
		dist = getDistance(data[cl2[k]], newCentroid2, d);
		if(dist < min) {
			min = dist;
			center2 = cl2[k];
		}
	}

	for(int j = 0; j < d; j++) {
		newCentroid1[j] = data[center1][j];
		newCentroid2[j] = data[center2][j];
	}
}

void KMedoidsSolver::swap(int p1, int p2) {
	int N = this->dataFrame.getInstance().N;
	int M = this->dataFrame.getInstance().M;
	int D = this->dataFrame.getInstance().D;
	int c1 = this->solution[p1];
	int c2 = this->solution[p2];
	double** data = this->dataFrame.getData();

	updateClusters(p1, p2, c1, c2);

	for(int j = 0; j < D; j++) {
		this->centroid[c1][j] = getMedian(this->clusters[c1], j);
		this->centroid[c2][j] = getMedian(this->clusters[c2], j);
	}

	double min;
	double dist;

	int centers[M];

	for(int i = 0; i < M; i++) {
		min = MAX_FLOAT;
		// setting the real centroid as the closest point to the artificial centroid calculated above
		for(int k = 0; k < this->clusters[i].size(); k++) {
			dist = getDistance(data[this->clusters[i][k]], this->centroid[i], D);
			if(dist < min) {
				min = dist;
				centers[i] = this->clusters[i][k];
			}
		}
	}

	for(int i = 0; i < M; i++) {
		for(int j = 0; j < D; j++) {
			this->centroid[i][j] = data[centers[i]][j];
		}
	}
}