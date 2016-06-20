#ifndef K_Medoids_Solver
#define K_Medoids_Solver

#include "KCenterSolver.h"

class KMedoidsSolver: public KCenterSolver {
	private:
		std::vector<int>* clusters;

	public:
		KMedoidsSolver(DataFrame dataFrame, int* solution);
		~KMedoidsSolver();
		std::vector<int>* getClusters() const;
		double getMedian(std::vector<int> cluster, int j);
		void calculateNewCentroids(int p, int c2, double* newCentroid1, double* newCentroid2);
		void calculateNewCentroidsSwap(int p1, int p2, double* newCentroid1, double* newCentroid2);
		void createCenters();
		void relocate(int p, int c2);
		void swap(int p1, int p2);
		void updateClusters(int p, int c1, int c2);
		void updateClusters(int p1, int p2, int c1, int c2);
};

#endif