/************************************************************************************
Solver.h
Solver

Created by Daniel Gribel

This header file contains the Solver class declaration.

The Solver class is an absract and super-class which is defined to solve the optimization (clustering) problem.
It is associated to a space of points defined in DataFrame and stores a solution, its cost, and
the cardinality of each cluster (group) of the solution.
*************************************************************************************/

#ifndef Solver_H
#define Solver_H

#include <algorithm>
#include <limits>
#include "Util.h"
#include "DataFrame.h"

class Solver {
	protected:
		
		/*Current solution the solver is working with*/
		int* solution;
		
		/*A reference to the data frame, which describes the dataset*/
		DataFrame dataFrame;
		
		/*The cost of the current solution*/
		double cost;

		/*The cardinality of each cluster of the solution*/
		int* cardinality;

	public:
		
		/*Get the current solution the solver is working with*/
		int* getSolution() { return this->solution; };
        
        /*Get the cost of the current solution*/
        double getCost() { return this->cost; };
        
		/*Get the cardinality of each cluster of the solution*/
        int* getCardinality() { return this->cardinality; };
        
        /*Get the data frame*/
        DataFrame getDataFrame() const { return this->dataFrame; };
        
        /*Performs the local search. This is one of the core parts of the programm, once it performs
        local improvements in the current solution. This method is implemented through polymorphism,
        according to how it is defined by classes that inherit Solver*/
        virtual void localSearch(std::vector<int>* conflictGraph) = 0;
		
		/*Calculates the solution cost. This method is implemented through polymorphism,
		according to how it is defined by classes that inherit Solver*/
		virtual void calculateCost() = 0;
};

#endif