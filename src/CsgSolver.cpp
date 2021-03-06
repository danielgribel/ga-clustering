/************************************************************************************
CsgSolver.cpp
CsgSolver

Created by Daniel Gribel

This cpp file contains the CsgSolver class definition.

The CsgSolver class represents an optimization solver that considers the total within
clusters distances as the objective for clustering.

Given an initial solution, it applies local improvements aiming to minimize the total
within group distances, i.e., the sum of distances of all points belonging to the
same cluster.
*************************************************************************************/

#include "CsgSolver.h"

using namespace std;

/*CsgSolver constructor*/
CsgSolver::CsgSolver(DataFrame* _dataFrame, int* _solution, string solverId) {
	setDataFrame(_dataFrame);
	setSolution(_solution);
	setSolverId(solverId);
	calculateCost();

	int n = this->dataFrame->getInstance().N;
	int m = this->dataFrame->getInstance().M;

	this->cardinality = new int[m];

	for(int i = 0; i < m; i++) {
		this->cardinality[i] = 0;
	}

	for(int i = 0; i < n; i++) {
		this->cardinality[solution[i]]++;
	}
}

/*CsgSolver destructor*/
CsgSolver::~CsgSolver() {

}

/*Get the distance of point p to every point within cluster c, i.e., the contribution of p on cluster c*/
double CsgSolver::costContribution(int p, int c) {
	int n = this->dataFrame->getInstance().N;
	double** sim = this->dataFrame->getSim();
	double cost = 0.0;

	for(int i = 0; i < n; i++) {
		if(this->solution[i] == c) {
			cost = cost + sim[i][p];
		}
	}

	return cost;
}

/*Performs the local search. This is one of the core parts of the programm, once it performs
local improvements in the current solution. In CsgSolver, the local search is perform in such a way
that moves aim to minimize the total within group distances*/
void CsgSolver::localSearch(std::vector<int>* conflictGraph) {
	double** data = this->dataFrame->getData();
	double** sim = this->dataFrame->getSim();
	int n = this->dataFrame->getInstance().N;
	int m = this->dataFrame->getInstance().M;
	bool improvingSolution = true;
	double newcost;
	int arr[n];
	int i;

	double increasing;
	double decreasing;

	int it = 0;
	bool sm;
	
	while(improvingSolution == true) {
		
		improvingSolution = false;

		/*Shuffle the order in which elements will be explored*/
		shuffle(arr, n); // O(n)

		for(int i1 = 0; i1 < n; i1++) {
			i = arr[i1];
			for(int k = 0; k < m; k++) {
				sm = shouldMove(conflictGraph[i], k, i);
				if((this->solution[i] != k) && (sm == true)) {

					/*Get the contribution of point p to its current cluster,
					as a value of decreasing in the total cost if the move is performed*/ 
					decreasing = costContribution(i, this->solution[i]);

					/*Get the contribution of point p to the destiny current cluster,
					as a value of increasing in the total cost if the move is performed*/
					increasing = costContribution(i, k);

					/*Calculate the solution new cost if the move is performed*/
					newcost = this->cost - decreasing + increasing;

					/*If the cost of the solution obtained by applying the move is less than the
					cost of the curent solution, then perform the move and update the current solution*/
					if(newcost < this->cost) {
						this->cardinality[solution[i]] = this->cardinality[solution[i]] - 1;
						this->cardinality[k] = this->cardinality[k] + 1;
						this->solution[i] = k;
						this->cost = newcost;
						improvingSolution = true;
					}
				}
			}
		}
		it++;
	}
}

/*Calculate solution cost from scratch. Heavy processing, should be called only once, when creating the Solver*/
void CsgSolver::calculateCost() {
	int n = this->dataFrame->getInstance().N;
	double** sim = this->dataFrame->getSim();
	this->cost = 0.0;

	for(int i = 0; i < n; i++) {
		for(int j = i+1; j < n; j++) {
			if(this->solution[i] == this->solution[j]) {
				this->cost = this->cost + sim[i][j];
			}
		}
	}
}

/*Verify the cost of a solution from scratch -- Useful for testing if the generated cost for the
best solution found is correct*/
double CsgSolver::verifyCost() {
	int n = this->dataFrame->getInstance().N;
	double** sim = this->dataFrame->getSim();
	double cst = 0.0;

	for(int i = 0; i < n; i++) {
		for(int j = i+1; j < n; j++) {
			if(this->solution[i] == this->solution[j]) {
				cst = cst + sim[i][j];
			}
		}
	}

	return cst;
}