#include "CsgSolver.h"

using namespace std;

CsgSolver::CsgSolver(DataFrame _dataFrame, int* _solution) {
	setDataFrame(_dataFrame);
	setSolution(_solution);
	calculateCost();
}

CsgSolver::~CsgSolver() {

}

void CsgSolver::setSolution(int* _solution) {
	this->solution = _solution;
}

void CsgSolver::setDataFrame(DataFrame _dataFrame) {
	this->dataFrame = _dataFrame;
}

bool CsgSolver::shouldMove(std::vector<int> conflicts, int destCluster, int p) {
	for(int q = 0; q < conflicts.size(); q++) {
		if(this->solution[conflicts[q]] == destCluster) {
			return false;
		}
	}

	// avoid that a cluster is left empty
	for(int i = 0; i < this->dataFrame.getInstance().N; i++) {
		if((this->solution[i] == this->solution[p]) && (i != p)) {
			return true;
		}
	}

	return false;
}

double CsgSolver::delta(int c, int p) {
	int n = this->dataFrame.getInstance().N;
	double** sim = this->dataFrame.getSim();
	double cost = 0.0;

	for(int i = 0; i < n; i++) {
		if(this->solution[i] == c) {
			cost = cost + sim[i][p];
		}
	}

	return cost;
}

void CsgSolver::localSearch(std::vector<int>* conflictGraph) {
	double** data = this->dataFrame.getData();
	double** sim = this->dataFrame.getSim();
	int n = this->dataFrame.getInstance().N;
	int m = this->dataFrame.getInstance().M;
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
		shuffle(arr, n); // O(n)

		for(int i1 = 0; i1 < n; i1++) {
			i = arr[i1];
			for(int k = 0; k < m; k++) {
				sm = shouldMove(conflictGraph[i], k, i);
				if((this->solution[i] != k) && (sm == true)) {

					decreasing = delta(this->solution[i], i);
					increasing = delta(k, i);
					newcost = this->cost - decreasing + increasing;

					if(newcost < this->cost) {
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

void CsgSolver::calculateCost() {
	int n = this->dataFrame.getInstance().N;
	double** sim = this->dataFrame.getSim();
	this->cost = 0.0;

	for(int i = 0; i < n; i++) {
		for(int j = i+1; j < n; j++) {
			if(this->solution[i] == this->solution[j]) {
				this->cost = this->cost + sim[i][j];
			}
		}
	}
}