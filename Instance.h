#ifndef INSTANCE
#define INSTANCE

#include <string>
#include <vector>

struct Instance {
	int N; // number of points (dataset size)
	int M; // number of clusters (classes)
	int D; // number of attributes
	int labelCol;
	char delimiter;
	std::string file;
	bool hasId;
};

Instance createInstance(int N, int M, int D, int labelCol, char delimiter, std::string file, bool hasId);
std::vector <Instance> getDatasets();

#endif