#include "DataFrame.h"

using namespace std;

DataFrame::DataFrame() {

}

DataFrame::DataFrame(double** data, double** sim, int* label, vector< vector<int> > closest, Instance instance) {
	this->data = data;
	this->sim = sim;
	this->label = label;
	this->closest = closest;
	this->instance = instance;
}

DataFrame::~DataFrame() {
	
}

double** DataFrame::getData() const {
	return this->data;
}

double** DataFrame::getSim() const {
	return this->sim;
}

int* DataFrame::getLabel() const {
	return this->label;
}

vector< vector<int> > DataFrame::getClosest() const {
	return this->closest;
}

Instance DataFrame::getInstance() const {
	return this->instance;
}