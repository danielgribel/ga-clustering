#include "KCenterSolver.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

const double MAX_FLOAT = std::numeric_limits<double>::max();

KCenterSolver::KCenterSolver(DataFrame dataFrame, int* solution) {
	setDataFrame(dataFrame);
	setSolution(solution);
}

KCenterSolver::~KCenterSolver() {

}

double** KCenterSolver::getCentroids() const {
	return this->centroid;
}

void KCenterSolver::setSolution(int* newSolution) {
	this->solution = newSolution;
}

void KCenterSolver::setDataFrame(DataFrame newDataFrame) {
	this->dataFrame = newDataFrame;
}

bool KCenterSolver::shouldMove(std::vector<int> conflicts, int destCluster, int p) {
	for(int q = 0; q < conflicts.size(); q++) {
		if(this->solution[conflicts[q]] == destCluster) {
			return false;
		}
	}

	/*avoid that a cluster is left empty*/
	if(this->cardinality[solution[p]] > 1) {
		return true;
	}

	return false;
}

void KCenterSolver::localSearch(std::vector<int>* conflictGraph) {
	double** data = this->dataFrame.getData();
	double** sim = this->dataFrame.getSim();
	int n = this->dataFrame.getInstance().N;
	int m = this->dataFrame.getInstance().M;
	int d = this->dataFrame.getInstance().D;
	bool improvingSolution = true;
	double newcost;
	int prev = 0;
	int arr[n];
	int i;
	int modifiedClusters[m];
	int modifiedClusters2[m];

	double newCentroid1[d];
	double newCentroid2[d];

	int clusterI;
	int clusterJ;

	for(int t = 0; t < m; t++) {
		modifiedClusters2[t] = 1;
	}
	int it = 0;
	bool sm;
	
	while(improvingSolution == true) {
		
		improvingSolution = false;

		for(int t = 0; t < m; t++) {
			modifiedClusters[t] = 0;
		}
		
		shuffle(arr, n); // O(n)

		for(int i1 = 0; i1 < n; i1++) {
			i = arr[i1];
			for(int k = 0; k < m; k++) {
				sm = shouldMove(conflictGraph[i], k, i);
				if((this->solution[i] != k) && (sm == true) &&
					((modifiedClusters2[this->solution[i]] == 1) || (modifiedClusters2[k] == 1))) { // O(n) --> KMeans
					
					prev = this->solution[i];
					updateCentroidsRelocate(i, k, newCentroid1, newCentroid2); // O(d) --> KMeans
					newcost = getRelocateCost(i, k, newCentroid1, newCentroid2); // O(nd) --> KMeans

					if(newcost < this->cost) {
						relocate(i, k); // O(d) --> KMeans
						this->solution[i] = k;
						this->cost = newcost;
						improvingSolution = true;
						modifiedClusters[prev] = 1;
						modifiedClusters[k] = 1;
					}
				}
			}

			vector<int> closest = this->getDataFrame().getClosest()[i];

			for(int j = 0; j < closest.size(); j++) {
				clusterI = this->solution[i];
				clusterJ = this->solution[closest[j]];

				if((clusterI != clusterJ) && ((modifiedClusters2[clusterI] == 1) || (modifiedClusters2[clusterJ] == 1))) {

					updateCentroidsSwap(i, closest[j], newCentroid1, newCentroid2);
					newcost = getSwapCost(i, closest[j], newCentroid1, newCentroid2);
						
					if(newcost < this->cost) {
						swap(i, closest[j]);
						this->solution[i] = clusterJ;
						this->solution[closest[j]] = clusterI;
						this->cost = newcost;
						improvingSolution = true;
						modifiedClusters[clusterI] = 1;
						modifiedClusters[clusterJ] = 1;
					}
				}
			}
		}

		for(int t = 0; t < m; t++) {
			modifiedClusters2[t] = modifiedClusters[t];
		}
		it++;
	}
}

double KCenterSolver::getRelocateCost(int p, int c2, double* newCentroid1, double* newCentroid2) {
	int n = this->dataFrame.getInstance().N;
	int d = this->dataFrame.getInstance().D;
	int c1 = this->solution[p];
	double** data = this->dataFrame.getData();

	double c = 0.0;

	double* oldCentroid1 = this->centroid[c1];
	double* oldCentroid2 = this->centroid[c2];

	this->centroid[c1] = newCentroid1;
	this->centroid[c2] = newCentroid2;

	this->solution[p] = c2;

	for(int i = 0; i < n; i++) {
		c = c + getDistance(data[i], this->centroid[this->solution[i]], d);	
	}

	this->centroid[c1] = oldCentroid1;
	this->centroid[c2] = oldCentroid2;

	this->solution[p] = c1;

	return c;
}

double KCenterSolver::getSwapCost(int p1, int p2, double* newCentroid1, double* newCentroid2) {
	int n = this->dataFrame.getInstance().N;
	int d = this->dataFrame.getInstance().D;
	int c1 = this->solution[p1];
	int c2 = this->solution[p2];
	double** data = this->dataFrame.getData();

	double c = 0.0;

	double* oldCentroid1 = this->centroid[c1];
	double* oldCentroid2 = this->centroid[c2];

	this->centroid[c1] = newCentroid1;
	this->centroid[c2] = newCentroid2;

	this->solution[p1] = c2;
	this->solution[p2] = c1;

	for(int i = 0; i < n; i++) {
		c = c + getDistance(data[i], this->centroid[this->solution[i]], d);	
	}

	this->centroid[c1] = oldCentroid1;
	this->centroid[c2] = oldCentroid2;

	this->solution[p1] = c1;
	this->solution[p2] = c2;

	return c;
}

void KCenterSolver::calculateCost() {
	int n = this->dataFrame.getInstance().N;
	int d = this->dataFrame.getInstance().D;
	double** data = this->dataFrame.getData();

	this->cost = 0.0;
	
	for(int i = 0; i < n; i++) {
		this->cost = this->cost + getDistance(data[i], this->centroid[solution[i]], d);
	}
}