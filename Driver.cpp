#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <map>
#include <ctime>
#include <cstdlib>
#include <set>
#include <queue>

#include "Solver.h"
#include "CsgSolver.h"
#include "KCenterSolver.h"
#include "KMeansSolver.h"
#include "KMediansSolver.h"
#include "KMedoidsSolver.h"
#include "MinAssignment.h"
#include "Heap.h"
#include "HashTable.h"
#include "LinkedList.h"

using namespace std;

/*Declaration of a pair of integers*/
#define pii pair<int, int>

/*Declaration of a pair of double and integer*/
#define pdi pair<double, int>

/*Declaration of a pair of double and pair of integer*/
#define pip pair<double, pii>

#define F first

#define S second

/*Declaration of the datasets path*/
#define INPUT_PATH "data/"

const double MAX_FLOAT = std::numeric_limits<double>::max();

const double PORC_ANNOTATED = 0.0;

int lastImp = 0;

Solver* createSolver(DataFrame dataFrame, int* solution);

/*Print a solution represented by a partitioning (array representation)*/
void printSolution(int* solution, int n) {
	for(int i = 0; i < n; i++) {
		cout << solution[i] << ",";
	}
	cout << endl;
}

double rand(int a, int b, int c, int d) {
	int total = a + b + c + d;
	double randIndex = 1.0*(a + d)/total;
	return randIndex;
}

double crand(int a, int b, int c, int d) {
	int total = a + b + c + d;
	double crandIndex = (a - (1.0*(b + a)*(c + a))/total)/((1.0*(b + a + c + a))/2 - (1.0*(b + a)*(c + a))/total);
	return crandIndex;
}

void evalSolution(int* solution, int* label, Instance instance, string algorithm) {
	int a = 0;
	int b = 0;
	int c = 0;

	/*Count the number of pairs that are in the same cluster under C and in the same class under C'*/
	for(int i = 0; i < instance.N; i++) {
		for(int j = i+1; j < instance.N; j++) {
			if( (solution[i] == solution[j]) && (label[i] == label[j]) ) {
				a++;
			}
		}
	}

	int counter_solution[instance.M];
	int counter_results[instance.M];
	map <int, int> mapLabel;

	for(int i = 0; i < instance.M; i++) {
		counter_solution[i] = 0;
		counter_results[i] = 0;
	}
	
	// this is inefficient. there are smarter ways to do it
	for(int i = 0; i < instance.N; i++) {
		mapLabel[label[i]] = 0;
	}

	for(int i = 0; i < instance.N; i++) {
		counter_solution[solution[i]]++;
		mapLabel[label[i]]++;	
	}

	for(int i = 0; i < instance.M; i++) {
		b = b + (counter_solution[i])*(counter_solution[i] - 1)/2;
	}

	for(map<int,int>::iterator it=mapLabel.begin(); it!=mapLabel.end(); ++it) {
    	c = c + (it->second)*(it->second - 1)/2;
	}

	b = b - a;
	c = c - a;
	int d = instance.N*(instance.N-1)/2.0 - a - b - c;

	cout << algorithm + " #crand: " << crand(a, b, c, d) << endl;
}

int posInVector(vector<string> labels, string s) {
	int pos = std::find(labels.begin(), labels.end(), s) - labels.begin();
	return pos;
}

/*Delete a matrix of dimension m (memory free)*/
void deleteMatrix(double** matrix, int m) {
	for(int i = 0; i < m; i++) {
		delete [] matrix[i];
	}
	delete [] matrix;
}

DataFrame load(Instance instance) {
	int n = instance.N;
	int m = instance.M;
	int d = instance.D;

	int* label = new int[n];
	double **data = new double*[n];
	double **sim = new double*[n];
	
	for(int i = 0; i < n; i++) {
		data[i] = new double[d];
		sim[i] = new double[n];
	}

	vector<string> labels;
	int pos;

	ifstream file( (INPUT_PATH + instance.file).c_str() );

	int firstColumn = 0;
	int lastColumn = d + 1;

	if(instance.hasId == true) {
		lastColumn++;
	}

    for(int row = 0; row < n; row++) {
        string line;
        getline(file, line);
        if ( !file.good() ) {
        	cout << "Error on file line reading" << endl;
            break;
        }

        stringstream iss(line);
        int j = 0;
        for (int col = firstColumn; col < lastColumn; col++) {
            
            string val;
            getline(iss, val, instance.delimiter);

            /*Skip if dataset has ids*/
        	if((instance.hasId == true) && (col==firstColumn)) {
        		continue;
        	}
        	
            if(col == instance.labelCol) {
            	pos = posInVector(labels, val);
				if(pos >= labels.size()) {
					labels.push_back(val);
				}
            	label[row] = pos;
            }

            /*Get attribute value*/
            else {
            	data[row][j] = atof(val.c_str());
            	j++;
            }
        }
    }

    vector<pdi> closest;
    set<int> x;

    /*Fill similarity matrix with the distances between points of dataset*/
    int k = 0;
    for(int i = 0; i < n; i++) {
    	for(int j = i+1; j < n; j++) {
    		double dist = getDistance(data[i], data[j], d);
    		sim[i][j] = dist;
    		sim[j][i] = dist;
			k++;
			closest.push_back(pdi(dist, i));
			closest.push_back(pdi(dist, j));
    	}
    }

    const int numClosest = 20;

    vector< vector<int> > closestObjects(n, vector<int>(numClosest));
    vector<int> v;

    for(int i = 0; i < n; i++) {
    	closestObjects[i] = kSmallestIndices(sim[i], n, numClosest);
    }

	DataFrame dataFrame(data, sim, label, closestObjects, instance);

	return dataFrame;
}

void imposeConflicts(int* annotated, int* label, vector<int>* conflictGraph, int m) {
	int a;
	int b;

	for(int i = 0; i < m; i++) {
		a = annotated[i];
		for(int j = 0; j < m; j++) {
			b = annotated[j];
			if(label[a] != label[b]) {
				conflictGraph[a].push_back(b);
			}
		}
	}
}

int * intdup(int const * src, size_t len) {
   int * p = (int *) malloc(len * sizeof(int));
   memcpy(p, src, len * sizeof(int));
   return p;
}

/*vector<int*> getPopulation(Instance instance, int popSize) {
	int n = instance.N;
	int m = instance.M;
	vector<int*> population(popSize);

	for(int i = 0; i < popSize; i++) {
		population[i] = getRandomSolution(n, m);
	}

	return population;
}*/

void crossoverX(int* p1, int* p2, int* offspring, int n, int m) {
	vector< vector<int> > groups(2*m);
	int markedCluster[2*m];
	int markedPoint[n];

	int keyOrder[m];

	for(int i = 0; i < m; i++) {
		keyOrder[i] = i*2;
	}

	for(int i = 0; i < n; i++) {
		groups[keyOrder[p1[i]]].push_back(i);
		groups[keyOrder[p2[i]]+1].push_back(i);
	}

	for(int i = 0; i < 2*m; i++) {
		markedCluster[i] = 0;
	}

	for(int i = 0; i < n; i++) {
		markedPoint[i] = 0;
		offspring[i] = -1;
	}

	int max;
	int maxIndex;
	int maxC = 0;
	int counter = 0;
	int p;
	int c = 0;
	bool b;

	while(counter < 2*m) {
		max = 0;
		maxIndex = 0;

		// getting the group with max cardinality
		for(int i = 0; i < 2*m; i++) {
			if((groups[i].size() > max) && (markedCluster[i] == 0)) {
				max = groups[i].size();
				maxIndex = i;
			}
		}
		markedCluster[maxIndex] = 1;

		b = false;

		for(int i = 0; i < groups[maxIndex].size(); i++) {
			p = groups[maxIndex][i];
			if(markedPoint[p] == 0) {
				offspring[p] = c;	
				markedPoint[p] = 1;
				b = true;
			}
		}

		if(b == true) {
			c++;
		}
		if(c > maxC) {
			maxC = c;
		}
		if(c == m) {
			c = 0;
		}
		counter++;
	}
	if(maxC < m) {
		//printSolution(offspring, n);
		for(int i = 0; i < n; i++) {
			offspring[i] = p1[i];
		}
	}
}

/*Given a solution represented by a partitioning (array), convert it to a centroids representation*/
double** solutionToCentroids(int* solution, const DataFrame dataFrame) {
    const int n = dataFrame.getInstance().N;
    const int m = dataFrame.getInstance().M;
    const int d = dataFrame.getInstance().D;
    double** data = dataFrame.getData();

    double** centroid = new double*[m];
    int* sizes = new int[m];

    for(int i = 0; i < m; i++) {
        centroid[i] = new double[d];
        for(int j = 0; j < d; j++) {
            centroid[i][j] = 0.0;
        }
        sizes[i] = 0;
    }

    for(int i = 0; i < n; i++) {
        sizes[solution[i]] = sizes[solution[i]]+1;
        for(int j = 0; j < d; j++) {
            centroid[solution[i]][j] = centroid[solution[i]][j] + data[i][j];
        }
    }
    
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < d; j++) {
            if(sizes[i] != 0) {
                centroid[i][j] = (1.0*centroid[i][j])/sizes[i];
            }
        }
    }
    delete [] sizes;

    return centroid;
}

/*Given a solution represented by its centroids, convert it to a partitioning (array) representation*/
void centroidsToSolution(int* solution, double** centroid, const DataFrame dataFrame) {    
    const int n = dataFrame.getInstance().N;
    const int m = dataFrame.getInstance().M;
    const int d = dataFrame.getInstance().D;
    double** data = dataFrame.getData();

    double mindist;
    double dist;

    for(int i = 0; i < n; i++) {
        mindist = MAX_FLOAT;
        for(int j = 0; j < m; j++) {
        	dist = getDistance(data[i], centroid[j], d);
            if(dist < mindist) {
                mindist = dist;
                solution[i] = j;
            }
        }
    }
}

double** assignment(double** c1, double** c2, int m, int d) {
    double** matrix = new double*[m];

    for(int i = 0; i < m; i++) {
        matrix[i] = new double[m];
    }

    for(int i = 0; i < m; i++) {
        for(int j = 0; j < m; j++) {
            matrix[i][j] = getDistance(c1[i], c2[j], d);
        }
    }

    return matrix;
}

bool allClustersPopulated(int* solution, int n, int m) {
    vector< vector<int> > v(m);
    for(int i = 0; i < n; i++) {
        v[solution[i]].push_back(i);
    }
    for(int i = 0; i < m; i++) {
        if(v[i].size() == 0) {
            return false;
        }
    }
    return true;
}

void crossover(int* offspring1, int* p1, int* p2, const DataFrame dataFrame) {
    const int n = dataFrame.getInstance().N;
    const int m = dataFrame.getInstance().M;
    const int d = dataFrame.getInstance().D;
    
    double** c1 = solutionToCentroids(p1, dataFrame);
    double** c2 = solutionToCentroids(p2, dataFrame);
    double** c3 = new double* [m];

    for(int i = 0; i < m; i++) {
        c3[i] = new double[d];
    }

    double** matrix = assignment(c1, c2, m, d);

    std::vector<long> matching = minAssignment(matrix, m);

    int r;

    bool validOffspring = false;

    while(!validOffspring) {
        for(int i = 0; i < m; i++) {
            r = rand() % 2;
            if(r == 0) {
                for(int j = 0; j < d; j++) {
                    c3[i][j] = c1[i][j];
                }
            } else {
                for(int j = 0; j < d; j++) {
                    c3[i][j] = c2[matching[i]][j];
                }
            }
        }
        centroidsToSolution(offspring1, c3, dataFrame);
        validOffspring = allClustersPopulated(offspring1, n, m);
    }

    deleteMatrix(c1, m);
    deleteMatrix(c2, m);
    deleteMatrix(c3, m);
    deleteMatrix(matrix, m);
}

void getCardinality(int* cardinality, int* solution, int n, int m) {
    for(int i = 0; i < m; i++) {
        cardinality[i] = 0;
    }
    for(int i = 0; i < n; i++) {
    	cardinality[solution[i]] = cardinality[solution[i]] + 1; 	
    }
}

vector<int*> getPopulation(const int popSize, DataFrame dataFrame) {
    std::vector<int*> population;
    const int n = dataFrame.getInstance().N;
    const int m = dataFrame.getInstance().M;

    for(int i = 0; i < popSize; i++) {
        int* p = new int[n];
        for(int q = 0; q < n; q++) {
            p[q] = rand() % m + 0;
        }
        population.push_back(p);
    }

    return population;
}


vector<int*> selectSurvivors(HeapPdi* costHeap,
	vector<int*> population,
	vector<double>& solutionCost,
	int sizePopulation,
    DataFrame dataFrame) {

    const int n = dataFrame.getInstance().N;
    const int m = dataFrame.getInstance().M;

    HashTable* table = new HashTable();
    int maxPopulation = population.size();
    vector<int*> newPopulation;
    int* discarded = new int[maxPopulation];
    int id;
    HeapPdi* heapInd = new HeapPdi();
    HeapPdi* heapClones = new HeapPdi(); 
    int* cardinality = new int[m];

    for(int i = 0; i < maxPopulation; i++) {
        getCardinality(cardinality, population[i], n, m);

        if(table->existItem(cardinality, solutionCost[i], m)) {
            heapClones->push_max(solutionCost[i], i);
        } else {
            Item * anItem = new Item;
            (*anItem).cost = solutionCost[i];
            (*anItem).cardinality = cardinality;
            (*anItem).next = NULL;
            table->insertItem(anItem, m);
            heapInd->push_max(solutionCost[i], i);
        }
        discarded[i] = 0;
    }
    delete [] cardinality;
    int j = 0;
    
    while((j < (maxPopulation-sizePopulation)) && (heapClones->getHeap().size() > 0)) {
        id = heapClones->front_max().second;
        heapClones->pop_max();
        discarded[id] = 1;
        j++;
    }

    while(j < (maxPopulation-sizePopulation)) {
        id = heapInd->front_max().second;
        heapInd->pop_max();
        discarded[id] = 1;
        j++;
    }

    HeapPdi* costHeapAux = new HeapPdi();
    int l = 0;

    for(int i = 0; i < maxPopulation; i++) {
        if(discarded[i] == 0) {
            newPopulation.push_back(population[i]);
            costHeapAux->push_min(solutionCost[i], l);
            l++;
        } else {
            delete [] population[i];
        }
    }
    costHeap->setHeap(costHeapAux->getHeap());
    solutionCost.resize(sizePopulation);

    for(int i = 0; i < costHeap->getHeap().size(); i++) {
        solutionCost[costHeap->getHeap()[i].second] = costHeap->getHeap()[i].first;
    }

    delete heapInd;
    delete heapClones;
    delete costHeapAux;
    delete [] discarded;
    delete table;

    return newPopulation;
}

vector<int*> diversifyPopulation(HeapPdi* costHeap,
    vector<int*> population,
    vector<double>& solutionCost,
    const int numKeep,
    const int numNew,
    unsigned short m,
    DataFrame dataFrame,
    vector<int>* conflictGraph) {

    solutionCost.resize(numKeep + numNew);    
    vector<int*> newPopulation;
    double cost;
    int id;

    HeapPdi* costHeapAux = new HeapPdi();

    for(int i = 0; i < numKeep; i++) {
        cost = costHeap->front_min().first;
        id = costHeap->front_min().second;
        costHeapAux->push_min(cost, i);
        costHeap->pop_min();
        newPopulation.push_back(population[id]);
        solutionCost[i] = cost;
    }

    vector<int*> randomIndividuals = getPopulation(numNew, dataFrame);

    for(int i = 0; i < randomIndividuals.size(); i++) {
    	Solver* ind = createSolver(dataFrame, randomIndividuals[i]);
		ind->localSearch(conflictGraph);
        newPopulation.push_back(ind->getSolution());
        cost = ind->getCost();
        costHeapAux->push_min(cost, numKeep+i);
        solutionCost[numKeep+i] = cost;
    }

    costHeap->setHeap(costHeapAux->getHeap());

    delete costHeapAux;

    return newPopulation;
}

/*int* tournamentSelection(vector<int*> pop, int k, int sizePopulation, DataFrame dataFrame) {
    int* bestSolution;
    double indCost;
    double bestCost = MAX_FLOAT;
    
    for(int i = 0; i < k; i++) {
    	Solver* ind = createSolver(dataFrame, pop[rand() % sizePopulation]);
        indCost = ind->getCost();
        if(indCost < bestCost) {
            bestSolution = ind->getSolution();
            bestCost = indCost;
        }
        delete ind;
    }
    
    return bestSolution;
}*/

int* tournamentSelection(vector<int*> pop, vector<double> solutionCost) {
    int* best = NULL;
    double indCost;
    double bestCost = MAX_FLOAT;
    int r;
    const int sizePopulation = pop.size();

    for(int i = 0; i < 2; i++) {
        r = rand() % sizePopulation;
        indCost = solutionCost[r];
        if(indCost < bestCost) {
            best = pop[r];
            bestCost = indCost;
        }
    }

    return best;
}

double verifyCost(int* solution, DataFrame dataFrame) {
	int N = dataFrame.getInstance().N;
	int M = dataFrame.getInstance().M;
	int D = dataFrame.getInstance().D;
	double** data = dataFrame.getData();

	double centroid[M][D];
	int sizes[M];

	for(int i = 0; i < M; i++) {
	    for(int j = 0; j < D; j++) {
	    	centroid[i][j] = 0.0;
	    }
	    sizes[i] = 0;
	}

	for(int i = 0; i < N; i++) {
		sizes[solution[i]] = sizes[solution[i]]+1;
		for(int j = 0; j < D; j++) {
			centroid[solution[i]][j] = centroid[solution[i]][j] + data[i][j];
		}
	}
	
	for(int i = 0; i < M; i++) {
		for(int j = 0; j < D; j++) {
			centroid[i][j] = (1.0*centroid[i][j])/sizes[i];	
		}
	}

	double cost = 0.0;

	for(int i = 0; i < N; i++) {
		cost = cost + getDistance(data[i], centroid[solution[i]], D);
	}

	return cost;
}

void demo(int seed, clock_t begin) {
	const Instance instance = getDatasets()[0];
	const int n = instance.N;
	const int m = instance.M;
	const int d = instance.D;
	const int NUM_ANNOTATED = PORC_ANNOTATED * n;
	
	srand(seed);

	int* annotated;
	vector<int> conflictGraph[n];

	DataFrame dataFrame = load(instance);
	annotated = randAnnotated(n, NUM_ANNOTATED);
	imposeConflicts(annotated, dataFrame.getLabel(), conflictGraph, NUM_ANNOTATED);

	const int sizePopulation = 20;
	const int maxPopulation = 200;
	const int itNoImprovement = 400;
	const int itDiv = 100;
	const int maxIt = 1000;

	int it = 0;
	int lastImprovement = 0;
	int lastDiv = 0;
	double bestCost = MAX_FLOAT;
	double delta = 0.000000000;
	
	vector<int*> population;
	HeapPdi* costHeap = new HeapPdi();
	double costS;
	int* bestSolution = new int[n];
	std::vector<double> solutionCost;

	population = getPopulation(20, dataFrame);
	
	for(int i = 0; i < population.size(); i++) {
		Solver* ind = createSolver(dataFrame, population[i]);
		ind->localSearch(conflictGraph);
		population[i] = ind->getSolution();
		costS = ind->getCost();
        costHeap->push_min(costS, i);
        solutionCost.push_back(costS);
		if(costS < bestCost) {
			bestCost = costS;
            std::copy(ind->getSolution(), ind->getSolution() + n, bestSolution);
		}
		printf("r(%d) = %.15f\n", i, costS);
	}

	while(((it-lastImprovement) < itNoImprovement) && (it < maxIt)) {
		
		int* offspring1 = new int[n];
		int* p1 = tournamentSelection(population, solutionCost);
		int* p2 = tournamentSelection(population, solutionCost);

		//crossoverX(p1, p2, offspring1, n, m);
		crossover(offspring1, p1, p2, dataFrame);

		Solver* off1 = createSolver(dataFrame, offspring1);
		off1->localSearch(conflictGraph);
		double off1Cost = off1->getCost();

		if(off1Cost < bestCost) {
			bestCost = off1Cost;
			std::copy(off1->getSolution(), off1->getSolution() + n, bestSolution);
			lastImprovement = it;
		}

		population.push_back(off1->getSolution());
		costHeap->push_min(off1Cost, population.size() - 1);
		solutionCost.push_back(off1Cost);

		if(population.size() > maxPopulation) {
            population = selectSurvivors(costHeap, population, solutionCost, sizePopulation, dataFrame);
        }

        if( ((it-lastImprovement) >= itDiv) && ((it-lastDiv) >= itDiv) ) {
            lastDiv = it;
            population = diversifyPopulation(costHeap, population, solutionCost, sizePopulation, 2*sizePopulation, m, dataFrame, conflictGraph);
            if(costHeap->front_min().first < bestCost) {
                lastImprovement = it;
                bestCost = costHeap->front_min().first;
                std::copy(population[costHeap->front_min().second], population[costHeap->front_min().second] + n, bestSolution);
            }
        }

		it++;
		printf("%d) %.15f %s %d %s %d \n", it, bestCost, "size_pop =", (int)population.size(), "last_imp =", lastImprovement);
	}

	double costVerified = verifyCost(bestSolution, dataFrame);
	printf("%.15f ", bestCost);
	//printSolution(bestSolution, n);
	printf("%.15f ", costVerified);
	
	delete [] annotated;
	delete [] bestSolution;
}

Solver* createSolver(DataFrame dataFrame, int* solution) {

    /*Create a Solver instance. The object to be created must be of a class that inherit Solver*/
	Solver* solver = new CsgSolver(dataFrame, solution);

	return solver;
}

int main() {
    
    /*Start to count the running time*/
	clock_t begin = clock();

	/*Call the main function, which will perform the genetic loop*/
    demo(1607, begin);
    
    /*Stop the running time clock*/
	double elapsedSecs = double(clock() - begin) / CLOCKS_PER_SEC;
	
    /*Print the total time of execution*/
    cout << elapsedSecs << endl;
	
    return 0;
}