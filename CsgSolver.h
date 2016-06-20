#ifndef Csg_Solver
#define Csg_Solver

#include "Solver.h"

class CsgSolver: public Solver {
	public:
		CsgSolver(DataFrame _dataFrame, int* _solution);
		~CsgSolver();
		void setSolution(int* _solution);
		void setDataFrame(DataFrame _dataFrame);
		bool shouldMove(std::vector<int> conflicts, int destCluster, int p);
		double delta(int c, int p);
		void localSearch(std::vector<int>* conflictGraph);
		void calculateCost();
};

#endif