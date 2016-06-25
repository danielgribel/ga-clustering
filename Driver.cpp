/************************************************************************************
Driver.h
Driver

Created by Daniel Gribel

The Driver file is the main file in the system. It contains the GA (Genetic Algorithm)
main loop and some important directives as input file loading, input validation,
constants declaration, accuracy calculation for results, etc
*************************************************************************************/

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

/*Stores the id of the chosen solver*/
string solverId;

/*Declaration of a pair of integers*/
#define pii pair<int, int>

/*Declaration of a pair of double and integer*/
#define pdi pair<double, int>

/*Declaration of a pair of double and pair of integer*/
#define pip pair<double, pii>

/*Declaration of the datasets path*/
#define INPUT_PATH "data/"

/*Solver id for KMeansSolver */
#define MEAN_ID "mean"

/*Solver id for KMediansSolver*/ 
#define MEDIAN_ID "median"

/*Solver id for KMedoidsSolver*/
#define MEDOID_ID "medoid"

/*Solver id for CsgSolver*/
#define CSG_ID "csg"

/*Max value for a float numer (used as infinite)*/
const double MAX_FLOAT = std::numeric_limits<double>::max();

/*The porcentage of items used as training set -- Used only for classification purpose*/
const double PORC_ANNOTATED = 0.0;

Solver* createSolver(DataFrame* dataFrame, int* solution);

/*Print a solution represented by a partitioning (array representation)*/
void printSolution(int* solution, int n) {
	for(int i = 0; i < n; i++) {
		cout << solution[i] << ",";
	}
	cout << endl;
}

/*Calculate the rand index, which is a measure for partitions agreement. Useful to test accuracy*/
double rand(int a, int b, int c, int d) {
	int total = a + b + c + d;
	double randIndex = 1.0*(a + d)/total;
	return randIndex;
}

/*Calculate the c-rand index (or adjusted rand), which is a measure for partitions agreement.
Useful to test accuracy*/
double crand(int a, int b, int c, int d) {
	int total = a + b + c + d;
	double crandIndex = (a - (1.0*(b + a)*(c + a))/total)/((1.0*(b + a + c + a))/2 - (1.0*(b + a)*(c + a))/total);
	return crandIndex;
}

void evalSolution(int* solution, int* label, DataFrame* dataFrame, string algorithm) {
	Instance instance = dataFrame->getInstance();

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

/*Returns the position where a string occurs in a vector*/
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

/*Create an instance for a dataset*/
Instance getInstance(int n, int m, int d, char adelimiter, string afile, bool ahasId) {
    Instance instance;
    instance.N = n;
    instance.M = m;
    instance.D = d;
    instance.delimiter = adelimiter;
    instance.file = afile;
    instance.hasId = ahasId;

    return instance;
}

/*Load the dataset*/
DataFrame* load(string fileName, int m) {
    
    ifstream file( (INPUT_PATH + fileName).c_str() );

    /*Abort if the file does not exist*/
    if ( !file.good() ) {
        cout << "Error on file reading" << endl;
        return NULL;
    }

    int n;
    int d;

    file >> n;
    file >> d;

    /*Abort if the number of points in the dataset is smaller than the number of desired clusters*/
    if(n < m) {
        cout << "n must be greater than m" << endl;
        return NULL;
    }

    const Instance instance = getInstance(n, m, d, ',', fileName, false);
    
    int* label = new int[n];
    double **data = new double*[n];
    double **sim = new double*[n];
    
    for(int i = 0; i < n; i++) {
        data[i] = new double[d];
        sim[i] = new double[n];
    }

    vector<string> labels;
    int pos;

    int firstColumn = 0;
    int lastColumn = d + 1;

    if(instance.hasId == true) {
        lastColumn++;
    }
    
    string line;
    getline(file, line);

    for(int row = 0; row < n; row++) {
        getline(file, line);

        if ( !file.good() ) {
            /*Abort if the file format is not ok*/
            cout << "Error on file line reading" << endl;
            return NULL;
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
            
            /*Get the label (class)*/
            if(col == d) {
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

    /*Get the $numCloset closest points for each point in the dataset*/
    for(int i = 0; i < n; i++) {
        closestObjects[i] = kSmallestIndices(sim[i], n, numClosest);
    }

    DataFrame* dataFrame = new DataFrame(data, sim, label, closestObjects, instance);
    
    return dataFrame;
}

/*Impose conflicts in the supervised learning version. This function keeps the points that must not be
in the same cluster, as we have labels (classes) for the training set -- Used only for classification purpose*/
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

/*Given a solution represented by a partitioning (array), convert it to a centroids representation*/
double** solutionToCentroids(int* solution, const DataFrame* dataFrame) {
    const int n = dataFrame->getInstance().N;
    const int m = dataFrame->getInstance().M;
    const int d = dataFrame->getInstance().D;
    double** data = dataFrame->getData();

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
void centroidsToSolution(int* solution, double** centroid, const DataFrame* dataFrame) {    
    const int n = dataFrame->getInstance().N;
    const int m = dataFrame->getInstance().M;
    const int d = dataFrame->getInstance().D;
    double** data = dataFrame->getData();

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

/*Build the assignment matrix for two solutions (each represented as a set of centroids), i.e.,
the distances between centroids of 2 different solutions */
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

/*Check if all clusters are populated, i.e., if no cluster is left empty.
If this function returns true, we have a valid solution.*/
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

/*Crossover: For each centroid of the offspring, we chose a random centroid among the centroids of the 2 parents*/
void crossover(int* offspring1, int* p1, int* p2, const DataFrame* dataFrame) {
    const int n = dataFrame->getInstance().N;
    const int m = dataFrame->getInstance().M;
    const int d = dataFrame->getInstance().D;
    
    double** c1 = solutionToCentroids(p1, dataFrame);
    double** c2 = solutionToCentroids(p2, dataFrame);
    double** c3 = new double* [m];

    for(int i = 0; i < m; i++) {
        c3[i] = new double[d];
    }

    double** matrix = assignment(c1, c2, m, d);

    /*Find the mininum assignment (full bipartite matching) -- here, dlib library function is called*/
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

/*Get the cardinality of each cluster in a given solution*/
void getCardinality(int* cardinality, int* solution, int n, int m) {
    for(int i = 0; i < m; i++) {
        cardinality[i] = 0;
    }
    for(int i = 0; i < n; i++) {
    	cardinality[solution[i]] = cardinality[solution[i]] + 1; 	
    }
}

/*Get a random population of size $popSize*/
vector<int*> getPopulation(const int popSize, DataFrame* dataFrame) {
    std::vector<int*> population;
    const int n = dataFrame->getInstance().N;
    const int m = dataFrame->getInstance().M;

    for(int i = 0; i < popSize; i++) {
        int* p = new int[n];
        for(int q = 0; q < n; q++) {
            p[q] = rand() % m + 0;
        }
        population.push_back(p);
    }

    return population;
}

/*Survivors selection: this function aims to select the best individuals to propagate. Its main purpose
is to keep the best individuals when the $maxPopulation is achieved. This procedure determines the $sizePopulation
individuals that will go on to the next generation, by discarding $maxPopulation - $sizePopulation individuals
that are either clones or bad regarding the fitness*/
vector<int*> selectSurvivors(HeapPdi* costHeap,
	vector<int*> population,
	vector<double>& solutionCost,
	int sizePopulation,
    DataFrame* dataFrame) {

    const int n = dataFrame->getInstance().N;
    const int m = dataFrame->getInstance().M;

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

        /*Check if solution already exist in population (clone detection)*/
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
    
    /*Remove clones while they exist*/
    while((j < (maxPopulation-sizePopulation)) && (heapClones->getHeap().size() > 0)) {
        id = heapClones->front_max().second;
        heapClones->pop_max();
        discarded[id] = 1;
        j++;
    }

    /*If $sizePopulation not achieved, remove individuals that are not clones but have bad fitness*/
    while(j < (maxPopulation-sizePopulation)) {
        id = heapInd->front_max().second;
        heapInd->pop_max();
        discarded[id] = 1;
        j++;
    }

    HeapPdi* costHeapAux = new HeapPdi();
    int l = 0;

    /*Get the survivors and delete the others*/
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

    /*Update solutions cost with remained individuals*/
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

/*Diversification: this function aims to ensure the diversity of the population. It is called whenever $itDiv
iterations happen without improving the best solution. It is performed by eliminating all but the best $numKeep
individuals of the population and creating 2x $numKeep new individuals as in the initialization phase
(randomly generated and then submitted to local search)*/
vector<int*> diversifyPopulation(HeapPdi* costHeap,
    vector<int*> population,
    vector<double>& solutionCost,
    const int numKeep,
    const int numNew,
    unsigned short m,
    DataFrame* dataFrame,
    vector<int>* conflictGraph) {

    solutionCost.resize(numKeep + numNew);    
    vector<int*> newPopulation;
    double cost;
    int id;

    HeapPdi* costHeapAux = new HeapPdi();

    /*Keep the $numKeep best individuals*/
    for(int i = 0; i < numKeep; i++) {
        cost = costHeap->front_min().first;
        id = costHeap->front_min().second;
        costHeapAux->push_min(cost, i);
        costHeap->pop_min();
        newPopulation.push_back(population[id]);
        solutionCost[i] = cost;
    }

    vector<int*> randomIndividuals = getPopulation(numNew, dataFrame);

    /*Generate $numNew new individuals (randomly)*/
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

/*Binary tournament selection: this function randomly selects 2 individuals (with uniform probability)
from the population and keeps the one with the best fitness to set as one of the parents. Then,
the same selection scheme is performed to set the second parent.*/
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

/*Verify the cost of a solution from scratch -- Useful for testing if the generated cost for the
best solution found is correct*/
double verifyCost(int* solution, DataFrame* dataFrame) {
	int N = dataFrame->getInstance().N;
	int M = dataFrame->getInstance().M;
	int D = dataFrame->getInstance().D;
	double** data = dataFrame->getData();

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

/*Check if the solver passed as argument if valid*/
bool validSolverId() {
    if(solverId != MEAN_ID
        && solverId != MEDIAN_ID
        && solverId != MEDOID_ID
        && solverId != CSG_ID) {
        
        cout << "Invalid solver. Please, enter a valid option: mean, median, medoid, csg." << endl;
        return false;
    }
    return true;
}

/*Check if the file passed as argument was successfully loaded*/
bool validDataFrame(DataFrame* dataFrame) {
    if(dataFrame == NULL) {
        return false;
    }
    return true;
}

/*Validate both solver and file passed as arguments*/
bool validInput(DataFrame* dataFrame) {
    if(validSolverId() && validDataFrame(dataFrame)) {
        return true;
    }
    return false;
}

/*The GA loop of the system. Initially, the method generates a random population of individuals.
Then, it applies successively a number of operators to select two parent individuals and combine them,
yielding a new individual (offspring), which is enhanced by a local search procedure and then added to
the population. Additionally, survivors selection and population diversification are applied when some
criteria are reached.*/
void demo(int seed, string fileName, int k) {
	   
    srand(seed);
    DataFrame* dataFrame = load(fileName, k);

    if(validInput(dataFrame)) {
        const int n = dataFrame->getInstance().N;
        const int m = dataFrame->getInstance().M;
        const int d = dataFrame->getInstance().D;
        const int NUM_ANNOTATED = PORC_ANNOTATED * n;

        int* annotated = randAnnotated(n, NUM_ANNOTATED);
        vector<int> conflictGraph[n];
        imposeConflicts(annotated, dataFrame->getLabel(), conflictGraph, NUM_ANNOTATED);

        /*The size of the population*/
        const int sizePopulation = 20;
        
        /*The maximum size a population can achieve (then survivors selection is performed)*/
        const int maxPopulation = 200;

        /*The number of iterations without improvements that the algorithm will run*/
        const int itNoImprovement = 400;

        /*The number of iterations without improvements determined to perform diversification*/
        const int itDiv = 100;

        /*The maximum number of iterations that the algorithm will run*/
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

        /*Generate the initial population -- individuals are randomly generated*/
        population = getPopulation(20, dataFrame);
        
        /*Submit initial population to education (local improvement)*/
        for(int i = 0; i < population.size(); i++) {
            Solver* ind = createSolver(dataFrame, population[i]);
            ind->localSearch(conflictGraph);
            population[i] = ind->getSolution();
            costS = ind->getCost();
            costHeap->push_min(costS, i);
            solutionCost.push_back(costS);

            /*Stores the solution if it is better than the best solution found so far*/
            if(costS < bestCost) {
                bestCost = costS;
                std::copy(ind->getSolution(), ind->getSolution() + n, bestSolution);
            }
            printf("r(%d) = %.15f\n", i, costS);
        }

        while(((it-lastImprovement) < itNoImprovement) && (it < maxIt)) {
            
            int* offspring1 = new int[n];

            /*Selects the first parent for crossover*/
            int* p1 = tournamentSelection(population, solutionCost);
            
            /*Selects the second parent for crossover*/
            int* p2 = tournamentSelection(population, solutionCost);

            /*Perform crossover: generate a offspring, given two parents p1 and p2*/
            crossover(offspring1, p1, p2, dataFrame);

            Solver* off1 = createSolver(dataFrame, offspring1);

            /*Offspring undergoes education (local improvement)*/
            off1->localSearch(conflictGraph);
            double off1Cost = off1->getCost();

            /*Stores offspring if it is better than the best solution found so far*/
            if(off1Cost < bestCost) {
                bestCost = off1Cost;
                std::copy(off1->getSolution(), off1->getSolution() + n, bestSolution);
                lastImprovement = it;
            }

            /*Add individual to population*/
            population.push_back(off1->getSolution());
            costHeap->push_min(off1Cost, population.size() - 1);
            solutionCost.push_back(off1Cost);

            /*If the size of population achieves $maxPopulation, then select survivors*/
            if(population.size() > maxPopulation) {
                population = selectSurvivors(costHeap, population, solutionCost, sizePopulation, dataFrame);
            }

            /*If $itDiv iterations happened without improving the best solution, then diversify population*/ 
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
        delete dataFrame;
    }
}

/*Create a Solver instance*/
Solver* createSolver(DataFrame* dataFrame, int* solution) {
    Solver* solver;

    /*Create a Solver instance. The object to be created must be of a class that inherit Solver*/
    if(solverId == MEAN_ID) {
        solver = new KMeansSolver(dataFrame, solution);
    } else if(solverId == MEDIAN_ID) {
        solver = new KMediansSolver(dataFrame, solution);
    } else if(solverId == MEDOID_ID) {
        solver = new KMedoidsSolver(dataFrame, solution);
    } else if(solverId == CSG_ID) {
        solver = new CsgSolver(dataFrame, solution);
    } else {
        solver = NULL;
        cout << "Invalid solver. Please, enter a valid option: mean, median, medoid, csg." << endl;
    }

	return solver;
}

int main(int argc, char** argv) {

    /*Get the arguments passed by the user*/
    string fileName = argv[1];
    solverId = argv[2];
    int m = atoi(argv[3]);
    
    //Start to count the running time
    clock_t begin = clock();

    //Call the main function, which will perform the genetic loop
    demo(1607, fileName, m);

    //Stop the running time clock
    double elapsedSecs = double(clock() - begin) / CLOCKS_PER_SEC;
        
    //Print the total time of execution
    cout << "Run time = " << elapsedSecs << endl;

    return 0;
}