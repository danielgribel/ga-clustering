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

using namespace std;

const double MAX_FLOAT = std::numeric_limits<double>::max();

/*KMedoidsSolver constructor*/
KMedoidsSolver::KMedoidsSolver(DataFrame* dataFrame, int* solution, string solverId)
: KCenterSolver(dataFrame, solution, solverId) {
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

/*Set the centroids. Given a solution, it calculates the centroids (medoid points within each cluster)*/
void KMedoidsSolver::createCenters() {
	int n = this->dataFrame->getInstance().N;
	int m = this->dataFrame->getInstance().M;
	int d = this->dataFrame->getInstance().D;

	double** data = this->dataFrame->getData();
	double** sim = this->dataFrame->getSim();

	this->centroid = new double*[m];
	this->clusters = new vector<int>[m];
	this->cardinality = new int[m];
	this->contrib = new double[m];

	double min;
	double dist;

	for(int i = 0; i < m; i++) {
		this->centroid[i] = new double[d];
		this->cardinality[i] = 0;
		this->contrib[i] = 0.0;
	}

	for(int i = 0; i < n; i++) {
		this->cardinality[solution[i]]++;
		this->clusters[solution[i]].push_back(i);
	}

	setDistances();

	int* bestMedoid = new int[m];

	for(int i = 0; i < m; i++) {
		min = MAX_FLOAT;
		for(int j = 0; j < this->clusters[i].size(); j++) {
			if(this->sumDist[i][j] < min) {
				bestMedoid[i] = this->clusters[i][j];
				min = this->sumDist[i][j];
			}
		}
		for(int k = 0; k < d; k++) {
			this->centroid[i][k] = data[bestMedoid[i]][k];
		}
	}

	int c;
	for(int i = 0; i < n; i++) {
		c = this->solution[i];
		this->contrib[c] = this->contrib[c] + sim[i][bestMedoid[c]];
	}

	delete [] bestMedoid;
}

/*//Get contrib of a point p in a new cluster c
double foo(int p, int c) {
	double contrib = 0.0;
	for(int i = 0; i < this->clusters[c].size(); i++) {
		contrib = contrib + sim[p][this->clusters[c][i]];
	}
}

double foo0(int p) {
	int c = this->solution[p];
	for(int i = 0; i < clusters[c].size(); i++) {
		if(p == clusters[c][i]) {
			return this->sumDist[c][i];
		}
	}
	return 0.0;
}*/

void KMedoidsSolver::setDistances() {
	int m = this->dataFrame->getInstance().M;
	int d = this->dataFrame->getInstance().D;
	double** sim = this->dataFrame->getSim();
	double sum = 0.0;

	this->sumDist = new vector<double>[m];

	for(int i = 0; i < m; i++) {
		for(int j = 0; j < this->clusters[i].size(); j++) {
			sum = 0.0;
			for(int k = 0; k < this->clusters[i].size(); k++) {
				sum = sum + sim[this->clusters[i][j]][this->clusters[i][k]];
			}
			this->sumDist[i].push_back(sum);
		}
	}
}

/*void KMedoidsSolver::setDistances(int c1, int c2) {
	double** sim = this->dataFrame->getSim();
	double sum = 0.0;

	for(int i = 0; i < this->clusters[c1].size(); i++) {
		sum = 0.0;
		for(int j = 0; j < this->clusters[c1].size(); j++) {
			sum = sum + sim[this->clusters[c1][i]][this->clusters[c1][j]];
		}
		this->sumDist[c1][i] = sum;
	}

	for(int i = 0; i < this->clusters[c2].size(); i++) {
		sum = 0.0;
		for(int j = 0; j < this->clusters[c2].size(); j++) {
			sum = sum + sim[this->clusters[c2][i]][this->clusters[c2][j]];
		}
		this->sumDist[c2][i] = sum;
	}
}*/

/*Update sumDist when performing relocate move*/
void KMedoidsSolver::setDistances(int q, int p, int c1, int c2) {
	double** sim = this->dataFrame->getSim();

	for(int i = 0; i < q; i++) {
		this->sumDist[c1][i] = this->sumDist[c1][i] - sim[this->clusters[c1][i]][p];
	}

	for(int i = q; i < this->sumDist[c1].size() - 1; i++) {
		this->sumDist[c1][i] = this->sumDist[c1][i+1] - sim[this->clusters[c1][i]][p];
	}

	int c1LastPosition = this->sumDist[c1].size() - 1;

	this->sumDist[c1].erase(this->sumDist[c1].begin() + c1LastPosition);

	double sum = 0.0;

	for(int i = 0; i < this->sumDist[c2].size(); i++) {
		this->sumDist[c2][i] = this->sumDist[c2][i] + sim[this->clusters[c2][i]][p];
		sum = sum + sim[this->clusters[c2][i]][p];
	}

	this->sumDist[c2].push_back(sum);
}

/*Update sumDist when performing swap move*/
void KMedoidsSolver::setDistances(int q1, int q2, int p1, int p2, int c1, int c2) {
	double** sim = this->dataFrame->getSim();
	this->sumDist[c1][q1] = 0.0;

	for(int i = 0; i < this->clusters[c1].size(); i++) {
		this->sumDist[c1][q1] = this->sumDist[c1][q1] + sim[p2][this->clusters[c1][i]];
	}

	for(int i = 0; i < this->clusters[c1].size(); i++) {
		if(i != q1) {
			this->sumDist[c1][i] = this->sumDist[c1][i] - sim[p1][this->clusters[c1][i]] + sim[p2][this->clusters[c1][i]];	
		}
	}

	// the same for c2
	this->sumDist[c2][q2] = 0.0;

	for(int i = 0; i < this->clusters[c2].size(); i++) {
		this->sumDist[c2][q2] = this->sumDist[c2][q2] + sim[p1][this->clusters[c2][i]];
	}

	for(int i = 0; i < this->clusters[c2].size(); i++) {
		if(i != q2) {
			this->sumDist[c2][i] = this->sumDist[c2][i] - sim[p2][this->clusters[c2][i]] + sim[p1][this->clusters[c2][i]];	
		}
	}

}

void KMedoidsSolver::updateContribs(int c1, int c2) {
	this->contrib[c1] = this->contrib1;
	this->contrib[c2] = this->contrib2;
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

	// here we should not discard
	//this->sumDist[c1].resize(this->clusters[c1].size());
	//this->sumDist[c2].resize(this->clusters[c2].size());

	updateContribs(c1, c2);
	setDistances(q, p, c1, c2);
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

	updateContribs(c1, c2);
	setDistances(q1, q2, p1, p2, c1, c2);
}

void KMedoidsSolver::updateCentroidsRelocate(int p, int c2, double* newCentroid1, double* newCentroid2) {
	int d = this->dataFrame->getInstance().D;
	int n = this->dataFrame->getInstance().N;
	int c1 = this->solution[p];

	double** data = this->dataFrame->getData();
	double** sim = this->dataFrame->getSim();
	
	int y = 0;

	for(int i = 0; i < this->clusters[c1].size(); i++) {
		if(this->clusters[c1][i] == p) {
			y = i;
		}
	}

	vector<int> cl1 = this->clusters[c1];
	vector<int> cl2 = this->clusters[c2];

	double min = MAX_FLOAT;

	for(int i = 0; i < sumDist[c1].size(); i++) {
		if(i != y && (sumDist[c1][i] - sim[cl1[i]][p] < min)) {
			this->bestMedoid1 = cl1[i];
			min = sumDist[c1][i] - sim[cl1[i]][p];
		}
	}

	min = 0.0;
	this->bestMedoid2 = p;
	
	//sum of distances from $p to every data point in c2
	for(int i = 0; i < sumDist[c2].size(); i++) {
		min = min + sim[p][cl2[i]];
	}

	for(int i = 0; i < sumDist[c2].size(); i++) {
		if(sumDist[c2][i] + sim[cl2[i]][p] < min) {
			this->bestMedoid2 = cl2[i];
			min = sumDist[c2][i] + sim[cl2[i]][p];
		}
	}

	for(int j = 0; j < d; j++) {
		newCentroid1[j] = data[this->bestMedoid1][j];
		newCentroid2[j] = data[this->bestMedoid2][j];
	}
}

void KMedoidsSolver::relocate(int p, int c2) {
	int d = this->dataFrame->getInstance().D;
	int c1 = this->solution[p];
	double** data = this->dataFrame->getData();
	double** sim = this->dataFrame->getSim();

	this->cardinality[c1] = this->cardinality[c1] - 1;
	this->cardinality[c2] = this->cardinality[c2] + 1;

	updateClusters(p, c1, c2);

	for(int j = 0; j < d; j++) {
		this->centroid[c1][j] = data[this->bestMedoid1][j];
		this->centroid[c2][j] = data[this->bestMedoid2][j];
	}
}

void KMedoidsSolver::updateCentroidsSwap(int p1, int p2, double* newCentroid1, double* newCentroid2) {
	int d = this->dataFrame->getInstance().D;
	int n = this->dataFrame->getInstance().N;
	int c1 = this->solution[p1];
	int c2 = this->solution[p2];

	double** data = this->dataFrame->getData();
	double** sim = this->dataFrame->getSim();

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

	double min = MAX_FLOAT;

	for(int i = 0; i < sumDist[c1].size(); i++) {
		if(i != y1 && (sumDist[c1][i] - sim[cl1[i]][p1] + sim[cl1[i]][p2] < min)) {
			this->bestMedoid1 = cl1[i];
			min = sumDist[c1][i] - sim[cl1[i]][p1] + sim[cl1[i]][p2];
		}
	}

	min = 0.0;
	this->bestMedoid2 = p1;
	
	//sum of distances from $p to every data point in c2
	for(int i = 0; i < sumDist[c2].size(); i++) {
		min = min + sim[p1][cl2[i]];
	}
	min = min - sim[p1][p2];

	for(int i = 0; i < sumDist[c2].size(); i++) {
		if(sumDist[c2][i] + sim[cl2[i]][p1] - sim[cl2[i]][p2] < min) {
			this->bestMedoid2 = cl2[i];
			min = sumDist[c2][i] + sim[cl2[i]][p1] - sim[cl2[i]][p2];
		}
	}

	for(int j = 0; j < d; j++) {
		newCentroid1[j] = data[this->bestMedoid1][j];
		newCentroid2[j] = data[this->bestMedoid2][j];
	}
}

/*Apply swap move between points p1 and p2 (p1 is moved to p2 cluster and p2 is moved to p1 cluster)*/
void KMedoidsSolver::swap(int p1, int p2) {
	int d = this->dataFrame->getInstance().D;
	int c1 = this->solution[p1];
	int c2 = this->solution[p2];
	double** data = this->dataFrame->getData();

	updateClusters(p1, p2, c1, c2);

	for(int j = 0; j < d; j++) {
		this->centroid[c1][j] = data[this->bestMedoid1][j];
		this->centroid[c2][j] = data[this->bestMedoid2][j];
	}
}

/*Get the solution cost if a relocate move is done*/
double KMedoidsSolver::getRelocateCost(int p, int c2, double* newCentroid1, double* newCentroid2) {
	int c1 = this->solution[p];

	double** sim = this->dataFrame->getSim();

	this->contrib1 = 0.0;
	this->contrib2 = 0.0;

	for(int i = 0; i < this->clusters[c1].size(); i++) {
		this->contrib1 = this->contrib1 + sim[bestMedoid1][clusters[c1][i]];
	}
	this->contrib1 = this->contrib1 - sim[bestMedoid1][p];

	for(int i = 0; i < this->clusters[c2].size(); i++) {
		this->contrib2 = this->contrib2 + sim[bestMedoid2][clusters[c2][i]];
	}
	this->contrib2 = this->contrib2 + sim[bestMedoid2][p];

	double c = this->cost - this->contrib[c1] - this->contrib[c2] + this->contrib1 + this->contrib2;

	return c;
}

/*Get the solution cost if a swap move is done*/
double KMedoidsSolver::getSwapCost(int p1, int p2, double* newCentroid1, double* newCentroid2) {
	int c1 = this->solution[p1];
	int c2 = this->solution[p2];

	double** sim = this->dataFrame->getSim();

	this->contrib1 = 0.0;
	this->contrib2 = 0.0;

	for(int i = 0; i < this->clusters[c1].size(); i++) {
		this->contrib1 = this->contrib1 + sim[bestMedoid1][clusters[c1][i]];
	}
	this->contrib1 = this->contrib1 - sim[bestMedoid1][p1] + sim[bestMedoid1][p2];

	for(int i = 0; i < this->clusters[c2].size(); i++) {
		this->contrib2 = this->contrib2 + sim[bestMedoid2][clusters[c2][i]];
	}
	this->contrib2 = this->contrib2 - sim[bestMedoid2][p2] + sim[bestMedoid2][p1];

	double c = this->cost - this->contrib[c1] - this->contrib[c2] + this->contrib1 + this->contrib2;

	return c;
}

/*Verify the cost of a solution from scratch -- Useful for testing if the generated cost for the
best solution found is correct*/
double KMedoidsSolver::verifyCost() {
	int n = this->dataFrame->getInstance().N;
	int m = this->dataFrame->getInstance().M;
	int d = this->dataFrame->getInstance().D;

	double** data = this->dataFrame->getData();
	double** sim = this->dataFrame->getSim();

	double** theCentroids = new double*[m];
	vector<int>* theClusters = new vector<int>[m];
	vector<double>* theSumDist = new vector<double>[m];
	int* theCardinality = new int[m];

	for(int i = 0; i < m; i++) {
		theCentroids[i] = new double[d];
		theCardinality[i] = 0;
	}

	for(int i = 0; i < n; i++) {
		theCardinality[this->solution[i]]++;
		theClusters[this->solution[i]].push_back(i);
	}

	double sum = 0.0;

	for(int i = 0; i < m; i++) {
		for(int j = 0; j < theClusters[i].size(); j++) {
			sum = 0.0;
			for(int k = 0; k < theClusters[i].size(); k++) {
				sum = sum + sim[theClusters[i][j]][theClusters[i][k]];
			}
			theSumDist[i].push_back(sum);
		}
	}

	int* bestMedoid = new int[m];
	double min;

	for(int i = 0; i < m; i++) {
		min = MAX_FLOAT;
		for(int j = 0; j < theClusters[i].size(); j++) {
			if(theSumDist[i][j] < min) {
				bestMedoid[i] = theClusters[i][j];
				min = theSumDist[i][j];
			}
		}
		for(int k = 0; k < d; k++) {
			theCentroids[i][k] = data[bestMedoid[i]][k];
		}
	}

	double theCost = 0.0;

	for(int i = 0; i < n; i++) {
		theCost = theCost + getDistance(data[i], data[bestMedoid[this->solution[i]]], d);
	}

	for(int i = 0; i < m; i++) {
		delete [] theCentroids[i];
	}
	delete [] theCentroids;
	delete [] theClusters;
	delete [] theSumDist;
	delete [] theCardinality;
	delete [] bestMedoid;

	return theCost;
}