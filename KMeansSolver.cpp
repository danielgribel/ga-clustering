#include "KMeansSolver.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

KMeansSolver::KMeansSolver(DataFrame dataFrame, int* solution)
: KCenterSolver(dataFrame, solution) {
	createCenters();
	calculateCost();
}

KMeansSolver::~KMeansSolver() {

}

double** KMeansSolver::getSumD() const {
	return this->sumD;
}

/*int* KMeansSolver::getSize() const {
	return this->size;
}*/

void KMeansSolver::calculateNewCentroids(int p, int c2, double* newCentroid1, double* newCentroid2) {
	double** data = this->dataFrame.getData();
	int c1 = this->solution[p];

	for(int i = 0; i < this->dataFrame.getInstance().D; i++) {
		newCentroid1[i] = (1.0*(this->sumD[c1][i] - data[p][i]))/(this->cardinality[c1]-1);
		newCentroid2[i] = (1.0*(this->sumD[c2][i] + data[p][i]))/(this->cardinality[c2]+1);
	}
}

void KMeansSolver::calculateNewCentroidsSwap(int p1, int p2, double* newCentroid1, double* newCentroid2) {
	double** data = this->dataFrame.getData();
	int c1 = this->solution[p1];
	int c2 = this->solution[p2];

	for(int i = 0; i < this->dataFrame.getInstance().D; i++) {
		newCentroid1[i] = (1.0*(this->sumD[c1][i] + data[p2][i] - data[p1][i]))/(this->cardinality[c1]);
		newCentroid2[i] = (1.0*(this->sumD[c2][i] + data[p1][i] - data[p2][i]))/(this->cardinality[c2]);
	}
}

void KMeansSolver::createCenters() {
	int N = this->dataFrame.getInstance().N;
	int M = this->dataFrame.getInstance().M;
	int D = this->dataFrame.getInstance().D;
	double** data = this->dataFrame.getData();
	
	this->centroid = new double*[M];
	this->sumD = new double*[M];
	this->cardinality = new int[M];

	for(int i = 0; i < M; i++) {
		this->centroid[i] = new double[D];
		this->sumD[i] = new double[D];
		this->cardinality[i] = 0;
	}

	for(int i = 0; i < M; i++) {
	    for(int j = 0; j < D; j++) {
	    	this->sumD[i][j] = 0.0;
	    }
	}

	for(int i = 0; i < N; i++) {
		this->cardinality[solution[i]]++;
		for(int j = 0; j < D; j++) {
			this->sumD[solution[i]][j] = this->sumD[solution[i]][j] + data[i][j];
		}
	}

	for(int i = 0; i < M; i++) {
		for(int j = 0; j < D; j++) {
			this->centroid[i][j] = (1.0*sumD[i][j])/cardinality[i];	
		}
	}
}

void KMeansSolver::relocate(int p, int c2) {
	int d = this->dataFrame.getInstance().D;
	int c1 = this->solution[p];
	double** data = this->dataFrame.getData();

	this->cardinality[c1] = this->cardinality[c1] - 1;
	this->cardinality[c2] = this->cardinality[c2] + 1;

	for(int i = 0; i < d; i++) {
		this->sumD[c1][i] = this->sumD[c1][i] - data[p][i];
		this->sumD[c2][i] = this->sumD[c2][i] + data[p][i];
		this->centroid[c1][i] = (1.0*this->sumD[c1][i])/this->cardinality[c1];
		this->centroid[c2][i] = (1.0*this->sumD[c2][i])/this->cardinality[c2];
	}
}

void KMeansSolver::swap(int p1, int p2) {
	int d = this->dataFrame.getInstance().D;
	int c1 = this->solution[p1];
	int c2 = this->solution[p2];
	double** data = this->dataFrame.getData();
	
	for(int i = 0; i < d; i++) {
		this->sumD[c1][i] = this->sumD[c1][i] + data[p2][i] - data[p1][i];
		this->sumD[c2][i] = this->sumD[c2][i] + data[p1][i] - data[p2][i];
		this->centroid[c1][i] = (1.0*this->sumD[c1][i])/this->cardinality[c1];
		this->centroid[c2][i] = (1.0*this->sumD[c2][i])/this->cardinality[c2];
	}
}