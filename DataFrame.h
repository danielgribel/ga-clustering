#ifndef DATA_FRAME
#define DATA_FRAME

#include <vector>

#include "Instance.h"

class DataFrame {
	private:
		double** data;
		double** sim;
		int* label;
		//std::vector<int> closest;
		std::vector< std::vector<int> > closest;
		Instance instance;

	public:
		DataFrame();
		//DataFrame(double** data, double** sim, int* label, std::vector<int> closest, Instance instance);
		DataFrame(double** data, double** sim, int* label, std::vector< std::vector<int> > closest, Instance instance);
		~DataFrame();
		double** getData() const;
		double** getSim() const;
		int* getLabel() const;
		//std::vector<int> getClosest() const;
		std::vector< std::vector<int> > getClosest() const;
		Instance getInstance() const;
};

#endif