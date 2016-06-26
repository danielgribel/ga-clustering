/************************************************************************************
CsgSolver.h
CsgSolver

Created by Daniel Gribel

This header file contains the CsgSolver class declaration.

The CsgSolver class represents an optimization solver that considers the total within
clusters distances as the objective for clustering.

Given an initial solution, it applies local improvements aiming to minimize the total
within group distances, i.e., the sum of distances of all points belonging to the
same cluster.
*************************************************************************************/

#ifndef Csg_Solver
#define Csg_Solver

#include "Solver.h"

class CsgSolver: public Solver {

	public:

		/*CsgSolver constructor*/
		CsgSolver(DataFrame* _dataFrame, int* _solution, std::string solverId);
		
		/*CsgSolver destructor*/
		~CsgSolver();

		/*Get the distance of point p to every point within cluster c, i.e., the contribution of p on cluster c*/
		double costContribution(int p, int c);
		
		/*Performs the local search. This is one of the core parts of the programm, once it performs
        local improvements in the current solution. In CsgSolver, the local search is perform in such a way
        that moves aim to minimize the total within group distances*/
		void localSearch(std::vector<int>* conflictGraph);
		
		/*Calculate solution cost from scratch. Heavy processing, should be called only once,
		when creating the Solver*/
		void calculateCost();

		/*Verify the cost of a solution from scratch -- Useful for testing if the generated cost for the
		best solution found is correct*/
		double verifyCost();
};

#endif