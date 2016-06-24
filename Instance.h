/************************************************************************************
Instance.h
Instance

Created by Daniel Gribel

This header file contains the Instance class declaration.
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

	/*In which column of the dataset file the label (class) is*/
	int labelCol;
	
	/*The delimiter for attributes (features) of the dataset (e.g.: comma, space, tab, etc)*/
	char delimiter;
	
	/*The file name of the dataset*/
	std::string file;

	/*Flag indicating if ids exist in the dataset. If True, the programm discard them, as they are not features*/
	bool hasId;
};

/*Create an instance for a dataset*/
Instance createInstance(int N, int M, int D, int labelCol, char delimiter, std::string file, bool hasId);

/*Add instances that will be loaded by the program and return them*/
std::vector <Instance> getDatasets();

#endif