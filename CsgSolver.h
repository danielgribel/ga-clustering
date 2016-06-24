/************************************************************************************
CsgSolver.h
CsgSolver

Created by Daniel Gribel

This header file contains the CsgSolver class declaration.
*************************************************************************************/

#ifndef Csg_Solver
#define Csg_Solver

#include "Solver.h"

class CsgSolver: public Solver {

	public:

		/*CsgSolver constructor*/
		CsgSolver(DataFrame _dataFrame, int* _solution);
		
		/*CsgSolver destructor*/
		~CsgSolver();

		/*Set a new solution*/
		void setSolution(int* _solution);

		/*Set a new data frame*/
		void setDataFrame(DataFrame _dataFrame);

		/*Check if is possible to perform a move. Possible reasons for move prohibition:
		- The move leaves a cluster empty
		- The move breaks some a-priori classification rule (when working with supervised classification)*/
		bool shouldMove(std::vector<int> conflicts, int destCluster, int p);
		
		/*Get the distance of point p to every point within cluster c, i.e., the contribution of p on cluster c*/
		double costContribution(int p, int c);
		
		/*Performs the local search. This is one of the core parts of the programm, once it performs
        local improvements in the current solution. In CsgSolver, the local search is perform in such a way
        that moves aim to minimize the total within group distances*/
		void localSearch(std::vector<int>* conflictGraph);
		
		/*Calculate solution cost from scratch. Heavy processing, should be called only once,
		when creating the Solver*/
		void calculateCost();
};

#endif