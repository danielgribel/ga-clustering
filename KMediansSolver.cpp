/************************************************************************************
KMediansSolver.cpp
KMediansSolver

Created by Daniel Gribel

This cpp file contains the KMediansSolver class definition.

The KMediansSolver class represents an optimization solver that considers the median point
of each cluster as a centroid. Given an initial solution, it applies local improvements
in order to minimize the distance of each point to the correspondent centroid -- in this
case, the median point inside the cluster.

The median point inside a cluster is obtained by calculating the median value for each
feature among all points belonging to the cluster. Thus, the centroid is not likely
to be a representative point (point belonging to the dataset).
*************************************************************************************/

#include "KMediansSolver.h"

using namespace std;

/*KMeansSolver constructor*/
KMediansSolver::KMediansSolver(DataFrame* dataFrame, int* solution)
: KCenterSolver(dataFrame, solution) {
	createCenters();
	calculateCost();
}

/*KMeansSolver destructor*/
KMediansSolver::~KMediansSolver() {

}

/*Get the list of elements belonging to each cluster*/
std::vector<int>* KMediansSolver::getClusters() const {
	return this->clusters;
}

/*Given the j-th cluster, get the median point, i.e., the median point for each feature,
which leads to a point that may not be a representative*/
double KMediansSolver::getMedian(std::vector<int> c, int j) {
	const int clusterSize = c.size();
	double values[clusterSize];
	double** data = this->dataFrame->getData();
	
	for(int i = 0; i < clusterSize; i++) {
		values[i] = data[c[i]][j];
	}

	int pos = clusterSize/2;

	nth_element(values, values + pos, values + clusterSize);
	double median = values[pos];
	
	return median;
}

/*Update clusters for relocate move. It removes point p from cluster c1 and add it to c2*/
void KMediansSolver::updateClusters(int p, int c1, int c2) {
	int q = 0;
	
	for(int i = 0; i < this->clusters[c1].size(); i++) {
		if(this->clusters[c1][i] == p) {
			q = i;
		}
	}

	this->clusters[c2].push_back(p);
	this->clusters[c1].erase(this->clusters[c1].begin() + q);
}

/*Update clusters for swap move. It adds point p1 to c2 and p2 to c1*/
void KMediansSolver::updateClusters(int p1, int p2, int c1, int c2) {
	int q1 = 0;
	int q2 = 0;
	
	for(int i = 0; i < this->clusters[c1].size(); i++) {
		if(this->clusters[c1][i] == p1) {
			q1 = i;
		}
	}

	for(int i = 0; i < this->clusters[c2].size(); i++) {
		if(this->clusters[c2][i] == p2) {
			q2 = i;
		}
	}

	this->clusters[c1][q1] = p2;
	this->clusters[c2][q2] = p1;
}

/*Set the centroids. Given a solution, it calculates the centroids (median points within each cluster)*/
void KMediansSolver::createCenters() {
	int N = this->dataFrame->getInstance().N;
	int M = this->dataFrame->getInstance().M;
	int D = this->dataFrame->getInstance().D;
	double** data = this->dataFrame->getData();

	this->centroid = new double*[M];
	this->clusters = new vector<int>[M];
	this->cardinality = new int[M];

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
}

/*Calculate the new centroids (median points within each cluster) if a relocate move is performed.
It does not change the current solution, but only obtain the new centroids that would be
resulted from a relocate move*/
void KMediansSolver::updateCentroidsRelocate(int p, int c2, double* newCentroid1, double* newCentroid2) {
	double** data = this->dataFrame->getData();
	int c1 = this->solution[p]; 
	int q = 0;
	int d = this->dataFrame->getInstance().D;

	for(int i = 0; i < this->clusters[c1].size(); i++) {
		if(this->clusters[c1][i] == p) {
			q = i;
		}
	}

	vector<int> cl1 = this->clusters[c1];
	vector<int> cl2 = this->clusters[c2];

	cl2.push_back(p);
	cl1.erase(cl1.begin() + q);

	for(int j = 0; j < d; j++) {
		newCentroid1[j] = getMedian(cl1, j);
		newCentroid2[j] = getMedian(cl2, j);
	}
}

/*Apply relocate move to point p (p is assigned to cluster c2)*/
void KMediansSolver::relocate(int p, int c2) {
	int N = this->dataFrame->getInstance().N;
	int M = this->dataFrame->getInstance().M;
	int D = this->dataFrame->getInstance().D;
	int c1 = this->solution[p];
	double** data = this->dataFrame->getData();

	this->cardinality[c1] = this->cardinality[c1] - 1;
	this->cardinality[c2] = this->cardinality[c2] + 1;

	updateClusters(p, c1, c2);

	for(int j = 0; j < D; j++) {
		this->centroid[c1][j] = getMedian(this->clusters[c1], j);
		this->centroid[c2][j] = getMedian(this->clusters[c2], j);
	}
}

/*Calculate the new centroids (median points within each cluster) if a swap move is performed.
It does not change the current solution, but only obtain the new centroids that would be
resulted from a swap move*/
void KMediansSolver::updateCentroidsSwap(int p1, int p2, double* newCentroid1, double* newCentroid2) {
	int D = this->dataFrame->getInstance().D;
	double** data = this->dataFrame->getData();
	int q1 = 0;
	int q2 = 0;
	int c1 = this->solution[p1];
	int c2 = this->solution[p2];

	for(int i = 0; i < this->clusters[c1].size(); i++) {
		if(this->clusters[c1][i] == p1) {
			q1 = i;
		}
	}

	for(int i = 0; i < this->clusters[c2].size(); i++) {
		if(this->clusters[c2][i] == p2) {
			q2 = i;
		}
	}

	vector<int> cl1 = this->clusters[c1];
	vector<int> cl2 = this->clusters[c2];

	cl1[q1] = p2;
	cl2[q2] = p1;

	for(int j = 0; j < D; j++) {
		newCentroid1[j] = getMedian(cl1, j);
		newCentroid2[j] = getMedian(cl2, j);
	}
}

/*Apply swap move between points p1 and p2 (p1 is moved to p2 cluster and p2 is moved to p1 cluster)*/
void KMediansSolver::swap(int p1, int p2) {
	int N = this->dataFrame->getInstance().N;
	int M = this->dataFrame->getInstance().M;
	int D = this->dataFrame->getInstance().D;
	int c1 = this->solution[p1];
	int c2 = this->solution[p2];
	double** data = this->dataFrame->getData();

	updateClusters(p1, p2, c1, c2);

	for(int j = 0; j < D; j++) {
		this->centroid[c1][j] = getMedian(this->clusters[c1], j);
		this->centroid[c2][j] = getMedian(this->clusters[c2], j);
	}
}