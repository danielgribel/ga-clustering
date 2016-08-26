/************************************************************************************
KCenterSolver.cpp
KCenterSolver

Created by Daniel Gribel

This cpp file contains the KCenterSolver class definition.

The KCenterSolver class represents a centroid-based optimization solver.
The KCenterSolver is an abstract class child of Solver and parent of
KMeansSolver, KMediansSolver and KMedoidsSolver classes, which have their own
implementation according to their objective functions.

In a general manner, given an initial solution, a centroid-based solver applies
local improvements in order to minimize the distance of each point to the
correspondent centroid -- KMeansSolver, KMediansSolver and KMedoidsSolver define how
the centroids are calculated.
*************************************************************************************/

#include "KCenterSolver.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

const double MAX_FLOAT = std::numeric_limits<double>::max();

/*KCenterSolver constructor*/
KCenterSolver::KCenterSolver(DataFrame* dataFrame, int* solution, std::string solverId) {
	setDataFrame(dataFrame);
	setSolution(solution);
	setSolverId(solverId);
}

/*KCenterSolver destructor*/
KCenterSolver::~KCenterSolver() {
	/*int m = this->dataFrame->getInstance().M;
	for(int i = 0; i < m; i++) {
		delete [] this->centroid[i];
	}
	delete [] this->centroid;*/
}

/*Get the list of centroids*/
double** KCenterSolver::getCentroids() const {
	return this->centroid;
}

/*Performs the local search. This is one of the core parts of the programm, once it performs
local improvements in the current solution. This method is implemented through polymorphism,
according to how it is defined by classes that inherit KCenterSolver*/
void KCenterSolver::localSearch(std::vector<int>* conflictGraph) {
	double** data = this->dataFrame->getData();
	double** sim = this->dataFrame->getSim();
	int n = this->dataFrame->getInstance().N;
	int m = this->dataFrame->getInstance().M;
	int d = this->dataFrame->getInstance().D;
	bool improvingSolution = true;
	double newcost;
	int prev = 0;
	int i;

	int* arr = new int[n];

	int* modifiedClusters = new int[m];
	int* modifiedClusters2 = new int[m];

	double* newCentroid1 = new double[d];
	double* newCentroid2 = new double[d];

	int clusterI;
	int clusterJ;

	for(int t = 0; t < m; t++) {
		modifiedClusters2[t] = 1;
	}
	int it = 0;
	bool sm;
	
	/*While any move improves the solution cost, keep doing it*/
	while(improvingSolution == true) {
		
		improvingSolution = false;

		for(int t = 0; t < m; t++) {
			modifiedClusters[t] = 0;
		}
		/*Shuffle the order in which elements will be explored*/
		shuffle(arr, n); // O(n)

		for(int i1 = 0; i1 < n; i1++) {
			i = arr[i1];
			for(int k = 0; k < m; k++) {
				sm = shouldMove(conflictGraph[i], k, i);

				/*Check if move is necessary*/
				if((this->solution[i] != k) && (sm == true) &&
					((modifiedClusters2[this->solution[i]] == 1) || (modifiedClusters2[k] == 1))) { // O(n) --> KMeans
					
					prev = this->solution[i];
					updateCentroidsRelocate(i, k, newCentroid1, newCentroid2); // O(d) --> KMeans
					newcost = getRelocateCost(i, k, newCentroid1, newCentroid2); // O(nd) --> KMeans

					/*If the cost of the solution obtained by applying the relocate move is less than the
					cost of the curent solution, then perform the relocate move and update the current solution*/
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

			vector<int> closest = this->getDataFrame()->getClosest()[i];

			for(int j = 0; j < closest.size(); j++) {
				clusterI = this->solution[i];
				clusterJ = this->solution[closest[j]];

				if((clusterI != clusterJ) && ((modifiedClusters2[clusterI] == 1) || (modifiedClusters2[clusterJ] == 1))) {

					updateCentroidsSwap(i, closest[j], newCentroid1, newCentroid2);
					newcost = getSwapCost(i, closest[j], newCentroid1, newCentroid2);
					
					/*If the cost of the solution obtained by applying the swap move is less than the
					cost of the curent solution, then perform the swap move and update the current solution*/ 	
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

	delete [] arr;
	delete [] modifiedClusters;
	delete [] modifiedClusters2;
	delete [] newCentroid1;
	delete [] newCentroid2;
}

/*Get the solution cost after a relocate move*/
double KCenterSolver::getRelocateCost(int p, int c2, double* newCentroid1, double* newCentroid2) {
	int n = this->dataFrame->getInstance().N;
	int d = this->dataFrame->getInstance().D;
	int c1 = this->solution[p];
	double** data = this->dataFrame->getData();

	double c = 0.0;

	/*Store current centroids*/
	double* oldCentroid1 = this->centroid[c1];
	double* oldCentroid2 = this->centroid[c2];

	/*Set new centroids*/
	this->centroid[c1] = newCentroid1;
	this->centroid[c2] = newCentroid2;

	/*Test relocation of point p to cluster c2*/
	this->solution[p] = c2;

	/*Calculate the new cost*/
	for(int i = 0; i < n; i++) {
		c = c + getDistance(data[i], this->centroid[this->solution[i]], d);	
	}

	/*Back to old centroids*/
	this->centroid[c1] = oldCentroid1;
	this->centroid[c2] = oldCentroid2;

	/*Send p back to old cluster*/
	this->solution[p] = c1;

	return c;
}

/*Get the solution cost after a swap move*/
double KCenterSolver::getSwapCost(int p1, int p2, double* newCentroid1, double* newCentroid2) {
	int n = this->dataFrame->getInstance().N;
	int d = this->dataFrame->getInstance().D;
	int c1 = this->solution[p1];
	int c2 = this->solution[p2];
	double** data = this->dataFrame->getData();

	double c = 0.0;

	/*Store current centroids*/
	double* oldCentroid1 = this->centroid[c1];
	double* oldCentroid2 = this->centroid[c2];

	/*Set new centroids*/
	this->centroid[c1] = newCentroid1;
	this->centroid[c2] = newCentroid2;

	/*Test swap of points*/
	this->solution[p1] = c2;
	this->solution[p2] = c1;

	/*Calculate the new cost*/
	for(int i = 0; i < n; i++) {
		c = c + getDistance(data[i], this->centroid[this->solution[i]], d);	
	}

	/*Back to old centroids*/
	this->centroid[c1] = oldCentroid1;
	this->centroid[c2] = oldCentroid2;

	/*Send points back to their old clusters*/
	this->solution[p1] = c1;
	this->solution[p2] = c2;

	return c;
}

/*Calculate solution cost from scratch*/
void KCenterSolver::calculateCost() {
	int n = this->dataFrame->getInstance().N;
	int d = this->dataFrame->getInstance().D;
	double** data = this->dataFrame->getData();

	this->cost = 0.0;

	for(int i = 0; i < n; i++) {
		this->cost = this->cost + getDistance(data[i], this->centroid[solution[i]], d);
	}
}