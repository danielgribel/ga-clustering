/************************************************************************************
KMeansSolver.h
KMeansSolver

Created by Daniel Gribel

This header file contains the KMeansSolver class declaration.
*************************************************************************************/

#ifndef K_Means_Solver
#define K_Means_Solver

#include "KCenterSolver.h"

class KMeansSolver: public KCenterSolver {

	private:

		/*The sum of each feature in each cluster. Useful to calculate the mean (average) point in the cluster*/
		double** sumFeatures;

	public:
		
		/*KMeansSolver constructor*/
		KMeansSolver(DataFrame data, int* solution);
		
		/*KMeansSolver destructor*/
		~KMeansSolver();
		
		/*Get the sum of each feature in each cluster */
		double** getSumFeatures() const;

		/*Calculate the new centroids (mean points within each cluster) if a relocate move is performed.
		It does not change the current solution, but only obtain the new centroids that would be
		resulted from a relocate move*/
		void updateCentroidsRelocate(int p, int c2, double* newCentroid1, double* newCentroid2);

		/*Calculate the new centroids (mean points within each cluster) if a swap move is performed.
		It does not change the current solution, but only obtain the new centroids that would be
		resulted from a swap move*/
		void updateCentroidsSwap(int p1, int p2, double* newCentroid1, double* newCentroid2);
		
		/*Set the centroids. Given a solution, it calculates the centroids (mean points within each cluster)*/
		void createCenters();

		/*Apply relocate move to point p (p is assigned to cluster c2)*/
		void relocate(int p, int c2);

		/*Apply swap move between points p1 and p2 (p1 is moved to p2 cluster and p2 is moved to p1 cluster)*/
		void swap(int p1, int p2);
};

#endif