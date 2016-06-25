/************************************************************************************
KMedoidsSolver.cpp
KMedoidsSolver

Created by Daniel Gribel

This cpp file contains the KMedoidsSolver class definition.

The KMedoidsSolver class represents an optimization solver that considers the medoid point
of each cluster as a centroid. Given an initial solution, it applies local improvements
in order to minimize the distance of each point to the correspondent centroid -- in this
case, the medoid point inside the cluster.

The medoid point inside a cluster is the representative (point belonging to the dataset)
that minimizes the distance for each other point within the cluster.
*************************************************************************************/

#include "KMedoidsSolver.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

const double MAX_FLOAT = std::numeric_limits<double>::max();

/*KMedoidsSolver constructor*/
KMedoidsSolver::KMedoidsSolver(DataFrame* dataFrame, int* solution)
: KCenterSolver(dataFrame, solution) {
	createCenters();
	calculateCost();
}

/*KMedoidsSolver destructor*/
KMedoidsSolver::~KMedoidsSolver() {

}

/*Get the list of elements belonging to each cluster*/
std::vector<int>* KMedoidsSolver::getClusters() const {
	return this->clusters;
}

/*Given the j-th cluster, get the median point, i.e., the median point for each feature,
which leads to a point that may not be a representative*/
double KMedoidsSolver::getMedian(std::vector<int> cl, int j) {
	const int clusterSize = cl.size();
	double values[clusterSize];
	double** data = this->dataFrame->getData();
	
	for(int i = 0; i < clusterSize; i++) {
		values[i] = data[cl[i]][j];
	}

	int pos = clusterSize/2;

	nth_element(values, values + pos, values + clusterSize);
	double median = values[pos];
	
	return median;
}

/*Update clusters for relocate move. It removes point p from cluster c1 and add it to c2*/
void KMedoidsSolver::updateClusters(int p, int c1, int c2) {
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
void KMedoidsSolver::updateClusters(int p1, int p2, int c1, int c2) {
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

/*Set the centroids. Given a solution, it calculates the centroids (medoid points within each cluster)*/
void KMedoidsSolver::createCenters() {
	int N = this->dataFrame->getInstance().N;
	int M = this->dataFrame->getInstance().M;
	int D = this->dataFrame->getInstance().D;
	double** data = this->dataFrame->getData();

	int* centers = new int[M];
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
		
		/*Setting the real centroid as the closest point to the artificial centroid calculated above*/
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

	delete [] centers;
}

/*Calculate the new centroids (medoid points within each cluster) if a relocate move is performed.
It does not change the current solution, but only obtain the new centroids that would be
resulted from a relocate move*/
void KMedoidsSolver::updateCentroidsRelocate(int p, int c2, double* newCentroid1, double* newCentroid2) {
	int d = this->dataFrame->getInstance().D;
	int n = this->dataFrame->getInstance().N;
	int c1 = this->solution[p];

	double** data = this->dataFrame->getData();
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

/*Apply relocate move to point p (p is assigned to cluster c2)*/
void KMedoidsSolver::relocate(int p, int c2) {
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

	double min;
	double dist;

	int centers[M];

	for(int i = 0; i < M; i++) {
		min = MAX_FLOAT;

		/*Setting the real centroid as the closest point to the artificial centroid calculated above*/
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

/*Calculate the new centroids (medoid points within each cluster) if a swap move is performed.
It does not change the current solution, but only obtain the new centroids that would be
resulted from a swap move*/
void KMedoidsSolver::updateCentroidsSwap(int p1, int p2, double* newCentroid1, double* newCentroid2) {
	int d = this->dataFrame->getInstance().D;
	int n = this->dataFrame->getInstance().N;
	int c1 = this->solution[p1];
	int c2 = this->solution[p2];

	double** data = this->dataFrame->getData();
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

/*Apply swap move between points p1 and p2 (p1 is moved to p2 cluster and p2 is moved to p1 cluster)*/
void KMedoidsSolver::swap(int p1, int p2) {
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

	double min;
	double dist;

	int centers[M];

	for(int i = 0; i < M; i++) {
		min = MAX_FLOAT;

		/*Setting the real centroid as the closest point to the artificial centroid calculated above*/
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