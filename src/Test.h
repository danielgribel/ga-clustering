/************************************************************************************
Test.h
Test

Created by Daniel Gribel

This header file contains the Test class declaration.

This class run some tests to check if the final solution generated is correct
*************************************************************************************/

#ifndef Test_H
#define Test_H

#include "Solver.h"

class Test {

	protected:
		
		/*Current solution*/
		Solver* solution;

		/*Small delta for cost precision*/
		double delta;
		
	public:

		/*Test constructor*/
		Test(Solver* asolution);
		
		/*Test destructor*/
		~Test();
		
		/*Verify the cost of a solution from scratch -- Useful for testing if the generated cost for the
		best solution found is correct*/
		bool verifyCost();

		/*Verify if any cluster is empty*/
		bool verifyCardinality();

		/*Run all tests*/
		void run();
};

#endif