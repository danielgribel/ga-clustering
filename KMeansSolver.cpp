/************************************************************************************
KMeansSolver.cpp
KMeansSolver

Created by Daniel Gribel

This cpp file contains the KMeansSolver class definition.
*************************************************************************************/

#include "KMeansSolver.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

/*KCenterSolver constructor*/
KMeansSolver::KMeansSolver(DataFrame dataFrame, int* solution)
: KCenterSolver(dataFrame, solution) {
	createCenters();
	calculateCost();
}

/*KCenterSolver destructor*/
KMeansSolver::~KMeansSolver() {

}

/*Get the sum of each feature in each cluster */
double** KMeansSolver::getSumFeatures() const {
	return this->sumFeatures;
}

/*Calculate the new centroids (mean points within each cluster) if a relocate move is performed.
It does not change the current solution, but only obtain the new centroids resulting from a relocate move*/
void KMeansSolver::updateCentroidsRelocate(int p, int c2, double* newCentroid1, double* newCentroid2) {
	double** data = this->dataFrame.getData();
	int c1 = this->solution[p];

	for(int i = 0; i < this->dataFrame.getInstance().D; i++) {
		newCentroid1[i] = (1.0*(this->sumFeatures[c1][i] - data[p][i]))/(this->cardinality[c1]-1);
		newCentroid2[i] = (1.0*(this->sumFeatures[c2][i] + data[p][i]))/(this->cardinality[c2]+1);
	}
}

/*Calculate the new centroids (mean points within each cluster) if a swap move is performed.
It does not change the current solution, but only obtain the new centroids resulting from a swap move*/
void KMeansSolver::updateCentroidsSwap(int p1, int p2, double* newCentroid1, double* newCentroid2) {
	double** data = this->dataFrame.getData();
	int c1 = this->solution[p1];
	int c2 = this->solution[p2];

	for(int i = 0; i < this->dataFrame.getInstance().D; i++) {
		newCentroid1[i] = (1.0*(this->sumFeatures[c1][i] + data[p2][i] - data[p1][i]))/(this->cardinality[c1]);
		newCentroid2[i] = (1.0*(this->sumFeatures[c2][i] + data[p1][i] - data[p2][i]))/(this->cardinality[c2]);
	}
}

/*Set the centroids. Given a solution, it calculates the centroids (mean points within each cluster)*/
void KMeansSolver::createCenters() {
	int N = this->dataFrame.getInstance().N;
	int M = this->dataFrame.getInstance().M;
	int D = this->dataFrame.getInstance().D;
	double** data = this->dataFrame.getData();
	
	this->centroid = new double*[M];
	this->sumFeatures = new double*[M];
	this->cardinality = new int[M];

	for(int i = 0; i < M; i++) {
		this->centroid[i] = new double[D];
		this->sumFeatures[i] = new double[D];
		this->cardinality[i] = 0;
	}

	for(int i = 0; i < M; i++) {
	    for(int j = 0; j < D; j++) {
	    	this->sumFeatures[i][j] = 0.0;
	    }
	}

	/*Sum features according to points belonging to the cluster*/
	for(int i = 0; i < N; i++) {
		this->cardinality[solution[i]]++;
		for(int j = 0; j < D; j++) {
			this->sumFeatures[solution[i]][j] = this->sumFeatures[solution[i]][j] + data[i][j];
		}
	}

	/*Create centroids based on points sum of features and cluster cardinality*/
	for(int i = 0; i < M; i++) {
		for(int j = 0; j < D; j++) {
			this->centroid[i][j] = (1.0*sumFeatures[i][j])/cardinality[i];	
		}
	}
}

/*Apply relocate move to point p (p is assigned to cluster c2)*/
void KMeansSolver::relocate(int p, int c2) {
	int d = this->dataFrame.getInstance().D;
	int c1 = this->solution[p];
	double** data = this->dataFrame.getData();

	/*Update clusters cardinality. c1 loses one point and c2 gains one.*/
	this->cardinality[c1] = this->cardinality[c1] - 1;
	this->cardinality[c2] = this->cardinality[c2] + 1;

	/*Update centroids according to the new set of points belonging to each cluster and their cardinalities*/ 
	for(int i = 0; i < d; i++) {
		this->sumFeatures[c1][i] = this->sumFeatures[c1][i] - data[p][i];
		this->sumFeatures[c2][i] = this->sumFeatures[c2][i] + data[p][i];
		this->centroid[c1][i] = (1.0*this->sumFeatures[c1][i])/this->cardinality[c1];
		this->centroid[c2][i] = (1.0*this->sumFeatures[c2][i])/this->cardinality[c2];
	}
}

/*Apply swap move between points p1 and p2 (p1 is moved to p2 cluster and p2 is moved to p1 cluster)*/
void KMeansSolver::swap(int p1, int p2) {
	int d = this->dataFrame.getInstance().D;
	int c1 = this->solution[p1];
	int c2 = this->solution[p2];
	double** data = this->dataFrame.getData();

	/*Update centroids according to the new set of points belonging to each cluster*/
	for(int i = 0; i < d; i++) {
		this->sumFeatures[c1][i] = this->sumFeatures[c1][i] + data[p2][i] - data[p1][i];
		this->sumFeatures[c2][i] = this->sumFeatures[c2][i] + data[p1][i] - data[p2][i];
		this->centroid[c1][i] = (1.0*this->sumFeatures[c1][i])/this->cardinality[c1];
		this->centroid[c2][i] = (1.0*this->sumFeatures[c2][i])/this->cardinality[c2];
	}
}