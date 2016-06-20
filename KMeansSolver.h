#ifndef K_Means_Solver
#define K_Means_Solver

#include "KCenterSolver.h"

class KMeansSolver: public KCenterSolver {
	private:
		double** sumD;
		//int* size;

	public:
		KMeansSolver(DataFrame data, int* solution);
		~KMeansSolver();
		double** getSumD() const;
		//int* getSize() const;
		void calculateNewCentroids(int p, int c2, double* newCentroid1, double* newCentroid2);
		void calculateNewCentroidsSwap(int p1, int p2, double* newCentroid1, double* newCentroid2);
		void createCenters();
		void relocate(int p, int c2);
		void swap(int p1, int p2);
};

#endif