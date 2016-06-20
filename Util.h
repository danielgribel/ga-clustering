#ifndef Util
#define Util

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>

struct greaterOp { bool operator()(std::pair<double,int> const& lhs, std::pair<double,int> const& rhs) { return lhs.first < rhs.first; } };
struct smallerOp { bool operator()(std::pair<double,int> const& lhs, std::pair<double,int> const& rhs) { return lhs.first > rhs.first; } };

void shuffle(int *myArray, size_t n);
int* getRandomSolution(int n, int m);
int* perturbation(int* s0, int k, int n, int m);
double manhattan(double* a, double* b, int d);
double euclidean(double* a, double* b, int d);
double getDistance(double* a, double* b, int d);
double* kLargest(double* myArray, int m, int k);
int* randAnnotated(int n, int m);
std::vector<int> kSmallestIndices(std::vector<double> v, int k);
std::vector<int> kSmallestIndices(double* v, int n, int k);

#endif