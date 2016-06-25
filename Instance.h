/************************************************************************************
Instance.h
Instance

Created by Daniel Gribel

This header file contains the Instance struct declaration.

Describes an instance (dataset) of the problem. It stores fields as number of points
(dataset size), number of desired clusters (groups), number of dimenions in the space
(features), the file name of the dataset, etc.
*************************************************************************************/

#ifndef INSTANCE
#define INSTANCE

#include <string>
#include <vector>

struct Instance {

	/*Number of points (dataset size)*/
	int N;

	/*Number of clusters (groups)*/
	int M;

	/*Number of dimenions in the space (features)*/
	int D;

	/*The delimiter for attributes (features) of the dataset (e.g.: comma, space, tab, etc)*/
	char delimiter;
	
	/*The file name of the dataset*/
	std::string file;

	/*Flag indicating if ids exist in the dataset. If True, the programm discard them, as they are not features*/
	bool hasId;
};

#endif