/************************************************************************************
Test.cpp
Test

Created by Daniel Gribel

This cpp file contains the Test class definition.

This class run some tests to check if the final solution generated is correct
*************************************************************************************/
#include <iostream>
#include <iomanip>

#include "Test.h"

using namespace std;

/*Test constructor*/
Test::Test(Solver* asolution) {
	this->solution = asolution;
	this->delta = 0.00001;
}

/*Verify the cost of a solution from scratch -- Useful for testing if the generated cost for the
best solution found is correct*/
bool Test::verifyCost() {

	if(this->solution->getCost() == this->solution->verifyCost()) {
		return true;
	}

	if((this->solution->getCost() > this->solution->verifyCost() - this->delta) &&
		(this->solution->getCost() < this->solution->verifyCost() + this->delta)) {
		return true;
	}

	cout << "this->solution->verifyCost() = " <<  setprecision(15) << this->solution->verifyCost() << endl;

	return false;
}

/*Verify if any cluster is empty*/
bool Test::verifyCardinality() {

	for(int i = 0; i < this->solution->getDataFrame()->getInstance().M; i++) {
		if(this->solution->getCardinality()[i] == 0) {
			return false;
		}
	}

	return true;
}

/*Run all tests*/
void Test::run() {

    DataFrame* dataFrame = this->solution->getDataFrame();

    cout << endl << ">> Tests:" << endl;

    cout << left << setw(dataFrame->getInstance().file.length()) << "Dataset" << "\t"
                << setw(6) << "n" << "\t"
                << setw(6) << "m" << "\t"
                << setw(10) << "Solver" << "\t"
                << setw(10) << "Cost" << "\t"
                << setw(10) << "Cardinality" << "\n";

    cout << left << setw(0) << dataFrame->getInstance().file << "\t"
                << setw(6) << dataFrame->getInstance().N << "\t"
                << setw(6) << dataFrame->getInstance().M << "\t"
                << setw(10) << this->solution->getSolverId() << "\t"
                << setw(10) << boolalpha << this->verifyCost() << "\t"
                << setw(10) << boolalpha << this->verifyCardinality() << "\n";

    cout << endl;
}