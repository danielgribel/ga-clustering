/************************************************************************************
DataFrame.cpp
DataFrame

Created by Daniel Gribel

This cpp file contains the DataFrame class definition.

Represents the data for an instance. It stores the dataset points in the space, the
similarity matrix between points of the dataset, the label (class) of points, etc.
*************************************************************************************/

#include "DataFrame.h"

using namespace std;

/*DataFrame constructor*/
DataFrame::DataFrame() {

}

/*DataFrame constructor*/
DataFrame::DataFrame(double** data, double** sim, int* label, vector< vector<int> > closest, Instance instance) {
	this->data = data;
	this->sim = sim;
	this->label = label;
	this->closest = closest;
	this->instance = instance;
}

/*DataFrame destructor*/
DataFrame::~DataFrame() {
	
}

/*Get the position of dataset points in the space*/
double** DataFrame::getData() const {
	return this->data;
}

/*Get the similarity matrix between points of the dataset*/
double** DataFrame::getSim() const {
	return this->sim;
}

/*Get label (class) of dataset points*/
int* DataFrame::getLabel() const {
	return this->label;
}

/*Get the list of closest points for each dataset point
This is required to apply swap moves in the local search phase*/
vector< vector<int> > DataFrame::getClosest() const {
	return this->closest;
}

/*Get the instance (dataset) being used*/
Instance DataFrame::getInstance() const {
	return this->instance;
}