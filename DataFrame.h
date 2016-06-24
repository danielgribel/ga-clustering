/************************************************************************************
DataFrame.h
DataFrame

Created by Daniel Gribel

This header file contains the DataFrame class declaration.
*************************************************************************************/

#ifndef DATA_FRAME
#define DATA_FRAME

#include <vector>

#include "Instance.h"

class DataFrame {
	
	private:

		/*Position of dataset points in the space: represented by a matrix nxd,
		where n is the number of points and d is the number of dimensions (features) in the space*/
		double** data;

		/*Similarity matrix between points of the dataset: represented by a nxn matrix,
		where n is the number of points*/
		double** sim;

		/*Label (class) of dataset points*/
		int* label;

		/*List of closest points for each dataset point
		Required to apply swap moves in the local search phase*/
		std::vector< std::vector<int> > closest;

		/*Instance (dataset) being used*/
		Instance instance;

	public:
		
		/*DataFrame constructor*/
		DataFrame();

		/*DataFrame constructor*/
		DataFrame(double** data, double** sim, int* label, std::vector< std::vector<int> > closest, Instance instance);
		
		/*DataFrame destructor*/
		~DataFrame();

		/*Get the position of dataset points in the space*/
		double** getData() const;

		/*Get the similarity matrix between points of the dataset*/
		double** getSim() const;

		/*Get label (class) of dataset points*/
		int* getLabel() const;

		/*Get the list of closest points for each dataset point
		This is required to apply swap moves in the local search phase*/
		std::vector< std::vector<int> > getClosest() const;

		/*Get the instance (dataset) being used*/
		Instance getInstance() const;
};

#endif