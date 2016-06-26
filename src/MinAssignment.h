/************************************************************************************
MinAssignment.h
MinAssignment

The contents of this file are in the public domain. See LICENSE_FOR_EXAMPLE_PROGRAMS.txt

This simple example shows how to call dlib's optimal linear assignment problem solver.
It is an implementation of the famous Hungarian algorithm and is quite fast, operating in
O(N^3) time.
*************************************************************************************/

#ifndef Min_Cost_Assignment
#define Min_Cost_Assignment

#include "../dlib-master/dlib/optimization/max_cost_assignment.h"
#include <iostream>

/*Returns the minimum assignment between two given sets S and S' with N objects each.
Used to find an assignment between centroids of two solutions*/
std::vector<long> minAssignment(double** mat, int m);

#endif