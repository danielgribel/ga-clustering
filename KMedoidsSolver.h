/************************************************************************************
KMedoidsSolver.h
KMedoidsSolver

Created by Daniel Gribel

This header file contains the KMedoidsSolver class declaration.

The KMedoidsSolver class represents an optimization solver that considers the medoid point
of each cluster as a centroid. Given an initial solution, it applies local improvements
in order to minimize the distance of each point to the correspondent centroid -- in this
case, the medoid point inside the cluster.

The medoid point inside a cluster is the representative (point belonging to the dataset)
that minimizes the distance for each other point within the cluster.
*************************************************************************************/

#ifndef K_Medoids_Solver
#define K_Medoids_Solver

#include "KCenterSolver.h"

class KMedoidsSolver: public KCenterSolver {
	
	private:

		/*Contains the list of elements belonging to each cluster*/
		std::vector<int>* clusters;

	public:

		/*KMedoidsSolver constructor*/
		KMedoidsSolver(DataFrame* dataFrame, int* solution);
		
		/*KMedoidsSolver destructor*/
		~KMedoidsSolver();

		/*Get the list of elements belonging to each cluster*/
		std::vector<int>* getClusters() const;

		/*Given the j-th cluster, get the median point, i.e., the median point for each feature,
		which leads to a point that may not be a representative*/
		double getMedian(std::vector<int> cluster, int j);

		/*Calculate the new centroids (medoid points within each cluster) if a relocate move is performed.
		It does not change the current solution, but only obtain the new centroids that would be
		resulted from a relocate move*/
		void updateCentroidsRelocate(int p, int c2, double* newCentroid1, double* newCentroid2);
		
		/*Calculate the new centroids (medoid points within each cluster) if a swap move is performed.
		It does not change the current solution, but only obtain the new centroids that would be
		resulted from a swap move*/
		void updateCentroidsSwap(int p1, int p2, double* newCentroid1, double* newCentroid2);

		/*Set the centroids. Given a solution, it calculates the centroids (medoid points within each cluster)*/
		void createCenters();

		/*Apply relocate move to point p (p is assigned to cluster c2)*/
		void relocate(int p, int c2);

		/*Apply swap move between points p1 and p2 (p1 is moved to p2 cluster and p2 is moved to p1 cluster)*/
		void swap(int p1, int p2);

		/*Update clusters for relocate move. It removes point p from cluster c1 and add it to c2*/
		void updateClusters(int p, int c1, int c2);
		
		/*Update clusters for swap move. It adds point p1 to c2 and p2 to c1*/
		void updateClusters(int p1, int p2, int c1, int c2);
};

#endif