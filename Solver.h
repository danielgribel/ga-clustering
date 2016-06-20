#ifndef Solver_H
#define Solver_H

#include <algorithm>
#include <limits>
#include "Util.h"
#include "DataFrame.h"

class Solver {
	protected:
		int* solution;
		DataFrame dataFrame;
		double cost;
		int* cardinality;

	public:
		int* getSolution() { return this->solution; };
        double getCost() { return this->cost; };
        int* getCardinality() { return this->cardinality; };
        DataFrame getDataFrame() const { return this->dataFrame; };
        virtual void localSearch(std::vector<int>* conflictGraph) = 0;
		virtual void calculateCost() = 0;
};

#endif