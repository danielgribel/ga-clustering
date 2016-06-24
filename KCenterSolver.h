#ifndef K_Center_Solver
#define K_Center_Solver

#include "Solver.h"

class KCenterSolver: public Solver {
	protected:
		double** centroid;

	public:
		KCenterSolver(DataFrame data, int* solution);
		~KCenterSolver();
		double** getCentroids() const;
		void setSolution(int* newSolution);
		void setDataFrame(DataFrame newDataFrame);
		bool shouldMove(std::vector<int> conflicts, int destCluster, int e);
		void localSearch(std::vector<int>* conflictGraph);
		virtual void updateCentroidsRelocate(int p, int c2, double* newCentroid1, double* newCentroid2) = 0;
		virtual void updateCentroidsSwap(int p1, int p2, double* newCentroid1, double* newCentroid2) = 0;
		virtual void createCenters() = 0;
		virtual void relocate(int p, int c2) = 0;
		virtual void swap(int p1, int p2) = 0;
		double getRelocateCost(int p, int c2, double* newCentroid1, double* newCentroid2);
		double getSwapCost(int p1, int p2, double* newCentroid1, double* newCentroid2);
		void calculateCost();
};

#endif