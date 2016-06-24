/************************************************************************************
Solver.cpp
Solver

Created by Daniel Gribel

This cpp file contains the Solver class definition.

The Solver class is an absract and super-class which is defined to solve the optimization (clustering) problem.
It is associated to a space of points defined in DataFrame and stores a solution, its cost, and
the cardinality of each cluster (group) of the solution.
*************************************************************************************/

#include "Solver.h"

/*Set a new solution*/
void Solver::setSolution(int* newSolution) {
	this->solution = newSolution;
}

/*Set a new data frame*/
void Solver::setDataFrame(DataFrame newDataFrame) {
	this->dataFrame = newDataFrame;
}

/*Check if is possible to perform a move. Possible reasons for move prohibition:
- The move leaves a cluster empty
- The move breaks some a-priori classification rule (when working with supervised classification)*/
bool Solver::shouldMove(std::vector<int> conflicts, int destCluster, int p) {
	for(int q = 0; q < conflicts.size(); q++) {
		if(this->solution[conflicts[q]] == destCluster) {
			return false;
		}
	}

	/*Avoid that a cluster is left empty*/
	if(this->cardinality[solution[p]] > 1) {
		return true;
	}

	return false;
}