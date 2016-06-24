/************************************************************************************
KCenterSolver.h
KCenterSolver

Created by Daniel Gribel

This header file contains the KCenterSolver class declaration.
*************************************************************************************/

#ifndef K_Center_Solver
#define K_Center_Solver

#include "Solver.h"

class KCenterSolver: public Solver {

	protected:

		/*The list of centroids, which are points in the space representing clusters.
		There is one centroid for each cluster*/
		double** centroid;

	public:

		/*KCenterSolver constructor*/
		KCenterSolver(DataFrame data, int* solution);
		
		/*KCenterSolver destructor*/
		~KCenterSolver();

		/*Get the list of centroids*/
		double** getCentroids() const;

		/*Set a new solution*/
		void setSolution(int* newSolution);
		
		/*Set a new data frame*/
		void setDataFrame(DataFrame newDataFrame);
		
		/*Check if is possible to perform a move. Possible reasons for move prohibition:
		- The move leaves a cluster empty
		- The move breaks some a-priori classification rule (when working with supervised classification)*/
		bool shouldMove(std::vector<int> conflicts, int destCluster, int e);
		
		/*Performs the local search. This is one of the core parts of the programm, once it performs
        local improvements in the current solution. This method is implemented through polymorphism,
        according to how it is defined by classes that inherit KCenterSolver*/
		void localSearch(std::vector<int>* conflictGraph);

		/*Calculate the new centroids if a relocate move is performed. It does not change the current solution,
		but only obtain the new centroids that would be resulted from a relocate move*/
		virtual void updateCentroidsRelocate(int p, int c2, double* newCentroid1, double* newCentroid2) = 0;
		
		/*Calculate the new centroids if a swap move is performed. It does not change the current solution,
		but only obtain the new centroids that would be resulted from a swap move*/
		virtual void updateCentroidsSwap(int p1, int p2, double* newCentroid1, double* newCentroid2) = 0;

		/*Set the centroids. This method is implemented through polymorphism, according to how
		centrality is defined by classes that inherit KCenterSolver*/
		virtual void createCenters() = 0;
		
		/*Apply relocate move to point p (p is assigned to cluster c2)*/
		virtual void relocate(int p, int c2) = 0;
		
		/*Apply swap move between points p1 and p2 (p1 is moved to p2 cluster and p2 is moved to p1 cluster)*/
		virtual void swap(int p1, int p2) = 0;

		/*Get the solution cost after a relocate move*/
		double getRelocateCost(int p, int c2, double* newCentroid1, double* newCentroid2);
		
		/*Get the solution cost after a swap move*/ 
		double getSwapCost(int p1, int p2, double* newCentroid1, double* newCentroid2);
		
		/*Calculate solution cost from scratch*/
		void calculateCost();
};

#endif