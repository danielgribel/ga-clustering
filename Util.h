/************************************************************************************
Util.h
Util

Created by Daniel Gribel

This header file contains the Util class declaration.
The Util class contains some useful functions for the system, as
distance calculations, k-largest indices, sequence randomization, etc
*************************************************************************************/

#ifndef Util
#define Util

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>

struct greaterOp { bool operator()(std::pair<double,int> const& lhs, std::pair<double,int> const& rhs) { return lhs.first < rhs.first; } };
struct smallerOp { bool operator()(std::pair<double,int> const& lhs, std::pair<double,int> const& rhs) { return lhs.first > rhs.first; } };

/*Shuffles an array: given an array of size n, returns a random permutation from 0 to n-1*/
void shuffle(int *myArray, size_t n);

/*Get a random solution with n elements and m clusters*/
int* getRandomSolution(int n, int m);

/*Apply perturbation to solution s0: set a random cluster (0..m) to k random elements*/
int* perturbation(int* s0, int k, int n, int m);

/*Calculate the Manhattan distance for two points a and b, given the number of dimensions d*/
double manhattan(double* a, double* b, int d);

/*Calculate the Euclidean distance for two points a and b, given the number of dimensions d*/
double euclidean(double* a, double* b, int d);

/*Calculate the distance between two points a and b, given the number of dimensions d
This function call the specific distance implementation (Manhattan, Euclidean, etc)*/
double getDistance(double* a, double* b, int d);

/*Given an array myArray, return the k largest elements*/
double* kLargest(double* myArray, int m, int k);

int* randAnnotated(int n, int m);

/*Given a vector v, return the indices of the k smallest elements*/
std::vector<int> kSmallestIndices(std::vector<double> v, int k);

/*Given an array v of size n, return the indices of the k smallest elements*/
std::vector<int> kSmallestIndices(double* v, int n, int k);

#endif