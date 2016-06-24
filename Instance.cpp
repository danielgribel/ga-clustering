/************************************************************************************
Instance.cpp
Instance

Created by Daniel Gribel

This cpp file contains the Instance class definition.
*************************************************************************************/

#include "Instance.h"

using namespace std;

/*Create an instance for a dataset*/
Instance createInstance(int N, int M, int D, int labelCol, char delimiter, string file, bool hasId) {
	Instance instance;
	instance.N = N;
	instance.M = M;
	instance.D = D;
	instance.labelCol = labelCol;
	instance.delimiter = delimiter;
	instance.file = file;
	instance.hasId = hasId;

	return instance;
}

/*Add instances that will be loaded by the program and return them*/
vector <Instance> getDatasets() {
	vector <Instance> instances;

	Instance iris45 = createInstance(45, 3, 4, 4, ',', "iris45.csv", false);
	Instance iris60 = createInstance(60, 3, 4, 4, ',', "iris60.csv", false);
	Instance iris = createInstance(150, 3, 4, 4, ',', "iris.csv", false);
	Instance fisher = createInstance(150, 3, 4, 4, ',', "fisher.csv", false);
	Instance bavaria1 = createInstance(89, 10, 3, 3, ',', "bavaria1.csv", false);
	Instance bavaria2 = createInstance(89, 10, 4, 4, ',', "bavaria2.csv", false);
	Instance wine = createInstance(178, 3, 13, 0, ',', "wine.csv", false);
	Instance glass = createInstance(214, 6, 9, 10, ',', "glass.csv", true);
	Instance heart = createInstance(297, 10, 13, 13, ',', "heart.csv", false);
	Instance liver = createInstance(345, 15, 6, 6, ',', "liver.csv", false);
	Instance ionosphere = createInstance(351, 25, 34, 34, ',', "ionosphere.csv", false);
	Instance breast = createInstance(699, 2, 9, 10, ',', "breast2.csv", true);
	Instance pima = createInstance(768, 10, 8, 8, ',', "pima.csv", false);
	Instance gesture = createInstance(1830, 5, 32, 32, ',', "a3_va3.csv", false);
	Instance spam = createInstance(4601, 2, 57, 57, ',', "spambase.csv", false);
	Instance yeast = createInstance(1484, 10, 8, 8, ',', "yeast.csv", false);

	instances.push_back(iris60);

	return instances;
}