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
#include "Test.h"

using namespace std;

/*Stores the id of the chosen solver*/
string solverId;

/*Flag to enable tests*/
#define RUN_TESTS true

/*Flag to enable accuracy calculation*/
bool accuracy = false;

/*Declaration of a pair of integers*/
#define pii pair<int, int>

/*Declaration of a pair of double and integer*/
#define pdi pair<double, int>

/*Declaration of a pair of double and pair of integer*/
#define pip pair<double, pii>

/*Declaration of the datasets path*/
#define INPUT_PATH "../data/"

/*Solver id for KMeansSolver */
#define MEAN_ID "mean"

/*Solver id for KMediansSolver*/ 
#define MEDIAN_ID "median"

/*Solver id for KMedoidsSolver*/
#define MEDOID_ID "medoid"

/*Solver id for CsgSolver*/
#define CSG_ID "csg"

#define EPS              (0.00000001)

/*Max value for a float numer (used as infinite)*/
const double MAX_FLOAT = std::numeric_limits<double>::max();

/*The porcentage of items used as training set -- Used only for classification purpose*/
const double PORC_ANNOTATED = 0.0;

Solver* createSolver(DataFrame* dataFrame, int* solution);

double evalSolution(int* solution, DataFrame* dataFrame);

/*Print a solution represented by a partition (array representation)*/
void printSolution(int* solution, int n) {
    cout << endl << ">> Best solution found:" << endl;
	for(int i = 0; i < n; i++) {
		cout << solution[i] << " ";
	}
    cout << endl;
}

/*Prints the summary of results*/
/*void printSummary(Solver* bestSolution, double elapsedSecs) {
    
    DataFrame* dataFrame = bestSolution->getDataFrame();
    
    printSolution(bestSolution->getSolution(), dataFrame->getInstance().N);

    cout << endl << ">> Summary of results:" << endl;

    cout << left << setw(dataFrame->getInstance().file.length()) << "Dataset" << "\t"
                << setw(2) << "n" << "\t"
                << setw(2) << "m" << "\t"
                << setw(8) << "Solver" << "\t"
                << setw(6) << "Obj function" << "\t"
                << setw(6) << "C-rand" << "\t"
                << setw(6) << "Time (s)" << "\n";

    cout << left << setw(0) << dataFrame->getInstance().file << "\t"
                << setw(2) << dataFrame->getInstance().N << "\t"
                << setw(2) << dataFrame->getInstance().M << "\t"
                << setw(8) << bestSolution->getSolverId() << "\t"
                << setw(6) << setprecision(8) << bestSolution->getCost() << "\t"
                << setw(6) << setprecision(4) << evalSolution(bestSolution->getSolution(), dataFrame) << "\t"
                << setw(6) << elapsedSecs << "\n";

    if(RUN_TESTS) {
        Test* test = new Test(bestSolution);
        test->run();
    }
}*/

void printSummary(Solver* bestSolution, double elapsedSecs) {
    
    DataFrame* dataFrame = bestSolution->getDataFrame();
    
    printSolution(bestSolution->getSolution(), dataFrame->getInstance().N);

    cout << endl << ">> Summary of results:" << endl;

    /*cout << left << dataFrame->getInstance().file << "\t"
                << dataFrame->getInstance().N << "\t"
                << dataFrame->getInstance().M << "\t"
                << bestSolution->getSolverId() << "\t"
                << bestSolution->getCost() << "\t"
                << evalSolution(bestSolution->getSolution(), dataFrame) << "\t"
                << elapsedSecs << " s\n";*/

    printf("%s\t", dataFrame->getInstance().file.c_str());
    printf("%d\t", dataFrame->getInstance().N);
    printf("%d\t", dataFrame->getInstance().M);
    printf("%s\t", bestSolution->getSolverId().c_str());
    printf("%.4f\t", bestSolution->getCost());
    printf("%.4f\t", evalSolution(bestSolution->getSolution(), dataFrame));
    printf("%.4f\n", elapsedSecs);

    if(RUN_TESTS) {
        Test* test = new Test(bestSolution);
        test->run();
        delete test;
    }
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

/*Evaluate solution accuracy*/
double evalSolution(int* solution, DataFrame* dataFrame) {
    if(!accuracy) {
        return 0;
    }

	Instance instance = dataFrame->getInstance();
    int* label = dataFrame->getLabel();

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

	return crand(a, b, c, d);
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
        cout << "Error: File does not exist or is corrupted." << endl;
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

    if(labels.size() == m) {
        accuracy = true;
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

/*Generate a random solution that assures that all clusters have been populated*/
void generateRandomPopulatedSolution(int* offspring1, const DataFrame* dataFrame) {
    const int n = dataFrame->getInstance().N;
    const int m = dataFrame->getInstance().M;
    const int d = dataFrame->getInstance().D;

    int* listClusters = new int[m];

    shuffle(listClusters, m);

    for(int i = 0; i < m; i++) {
        offspring1[i] = listClusters[i];
    }

    for(int i = m; i < n; i++) {
        offspring1[i] = rand() % m;
    }

    delete [] listClusters;
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
    int cont = 0;

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
        //std::cout << "validOffspring = " << validOffspring << std::endl;
        cont++;
        if(cont > 10) {
            generateRandomPopulatedSolution(offspring1, dataFrame);
            validOffspring = allClustersPopulated(offspring1, n, m);
        }
    }

    deleteMatrix(c1, m);
    deleteMatrix(c2, m);
    deleteMatrix(c3, m);
    deleteMatrix(matrix, m);
}

/*Swap (exchange) the assignment for two data points at positions i and j*/
void swap(int* v, const int i, const int j) {
    int t;
    t = v[i];
    v[i] = v[j];
    v[j] = t;
}

/*Mutate a solution with a rounding swap*/
void mutation(int* s0, int* s, int k, int n) {
    int* randItems = new int[n];
    for(int i = 0; i < n; i++) {
        s[i] = s0[i];
    }
    shuffle(randItems, n);
    for(int i = 0; i < 2*k; i = i+2) {
        swap(s, randItems[i], randItems[i+1]);
    }
    delete [] randItems;
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
vector<Solver*> getPopulation(const int popSize, DataFrame* dataFrame) {
    std::vector<Solver*> population;
    const int n = dataFrame->getInstance().N;
    const int m = dataFrame->getInstance().M;

    for(int i = 0; i < popSize; i++) {
        int* assignment = new int[n];
        for(int j = 0; j < n; j++) {
            assignment[j] = rand() % m + 0;
        }
        Solver* solution = createSolver(dataFrame, assignment);
        population.push_back(solution);
    }

    return population;
}

/*Survivors selection: this function aims to select the best individuals to propagate. Its main purpose
is to keep the best individuals when the $maxPopulation is achieved. This procedure determines the $sizePopulation
individuals that will go on to the next generation, by discarding $maxPopulation - $sizePopulation individuals
that are either clones or bad regarding the fitness*/
vector<Solver*> selectSurvivors(HeapPdi* costHeap,
	vector<Solver*> population,
	int sizePopulation,
    DataFrame* dataFrame) {

    const int n = dataFrame->getInstance().N;
    const int m = dataFrame->getInstance().M;

    HashTable* table = new HashTable();
    int maxPopulation = population.size();
    vector<Solver*> newPopulation;
    int* discarded = new int[maxPopulation];
    int id;
    HeapPdi* heapInd = new HeapPdi();
    HeapPdi* heapClones = new HeapPdi(); 

    //Add the best solution
    int topId = costHeap->front_min().second;
    double topCost = costHeap->front_min().first;
    Item * anItem = new Item;
    (*anItem).cost = topCost;
    (*anItem).cardinality = population[topId]->getCardinality();
    (*anItem).next = NULL;
    table->insertItem(anItem, m);
    heapInd->push_max(topCost, topId);

    for(int i = 0; i < maxPopulation; i++) {
        if(i != topId) {
            int* cardinality = population[i]->getCardinality();
            /*Check if solution already exist in population (clone detection)*/
            if(table->existItem(cardinality, population[i]->getCost(), m)) {
                heapClones->push_max(population[i]->getCost(), i);
            } else {
                Item * anItem = new Item;
                (*anItem).cost = population[i]->getCost();
                (*anItem).cardinality = cardinality;
                (*anItem).next = NULL;
                table->insertItem(anItem, m);
                heapInd->push_max(population[i]->getCost(), i);
            }    
        }
        discarded[i] = 0;
    }
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
            costHeapAux->push_min(population[i]->getCost(), l);
            l++;
        } else {
            delete population[i];
        }
    }
    costHeap->setHeap(costHeapAux->getHeap());

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
vector<Solver*> diversifyPopulation(HeapPdi* costHeap,
    vector<Solver*> population,
    const int numKeep,
    const int numNew,
    unsigned short m,
    DataFrame* dataFrame,
    vector<int>* conflictGraph) {

    vector<Solver*> newPopulation;
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
    }

    vector<Solver*> randomIndividuals = getPopulation(numNew, dataFrame);

    /*Generate $numNew new individuals (randomly)*/
    for(int i = 0; i < randomIndividuals.size(); i++) {
        randomIndividuals[i]->localSearch(conflictGraph);
        newPopulation.push_back(randomIndividuals[i]);
        cost = randomIndividuals[i]->getCost();
        costHeapAux->push_min(cost, numKeep+i);
    }

    costHeap->setHeap(costHeapAux->getHeap());

    delete costHeapAux;

    return newPopulation;
}

/*Binary tournament selection: this function randomly selects 2 individuals (with uniform probability)
from the population and keeps the one with the best fitness to set as one of the parents. Then,
the same selection scheme is performed to set the second parent.*/
int* tournamentSelection(vector<Solver*> pop) {
    Solver* best = NULL;
    double indCost;
    double bestCost = MAX_FLOAT;
    int r;
    const int sizePopulation = pop.size();

    for(int i = 0; i < 2; i++) {
        r = rand() % sizePopulation;
        indCost = pop[r]->getCost();
        if(indCost < bestCost) {
            best = pop[r];
            bestCost = indCost;
        }
    }

    return best->getSolution();
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
    
    //Start to count the running time
    clock_t begin = clock();
    
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
        const int itNoImprovement = 500;

        /*The number of iterations without improvements determined to perform diversification*/
        const int itDiv = 200;

        /*The number of exchanges on rounding mutation*/
        const int mutationStrength = 20;

        /*The maximum number of iterations that the algorithm will run*/
        const int maxIt = 1000;

        int it = 0;
        int lastImprovement = 0;
        int lastDiv = 0;
        double delta = 0.000000000;
        
        vector<Solver*> population;
        HeapPdi* costHeap = new HeapPdi();
        double costS;

        /*Generate the initial population -- individuals are randomly generated*/
        population = getPopulation(20, dataFrame);
        Solver* bestSolution = createSolver(dataFrame, population[0]->getSolution());

        /*Submit initial population to education (local improvement)*/
        for(int i = 0; i < population.size(); i++) {
            population[i]->localSearch(conflictGraph);
            costHeap->push_min(population[i]->getCost(), i);

            /*Stores the solution if it is better than the best solution found so far*/
            if(population[i]->getCost() < bestSolution->getCost()) {
                bestSolution = population[i];
            }

            printf("%d) = %.15f\n", i, population[i]->getCost());
        }

        while(it < 1000) {
        //while(((it-lastImprovement) < itNoImprovement) && (it < maxIt)) {
            
            int* offspring1 = new int[n];
            int* offspring2 = new int[n];

            /*Selects the first parent for crossover*/
            int* p1 = tournamentSelection(population);
            
            /*Selects the second parent for crossover*/
            int* p2 = tournamentSelection(population);

            /*Perform crossover: generate a offspring, given two parents p1 and p2*/
            crossover(offspring1, p1, p2, dataFrame);

            Solver* off1 = createSolver(dataFrame, offspring1);

            /*Offspring undergoes education (local improvement)*/
            off1->localSearch(conflictGraph);

            /*Stores offspring if it is better than the best solution found so far*/
            if(off1->getCost() < bestSolution->getCost()) {
                bestSolution = off1;
                lastImprovement = it;
            }

            /*Apply mutation on offspring*/
            mutation(off1->getSolution(), offspring2, mutationStrength, n);

            Solver* off2 = createSolver(dataFrame, offspring2);
            //off2->localSearch(conflictGraph);

            if(off2->getCost() < bestSolution->getCost()) {
                bestSolution = off2;
                lastImprovement = it;
            }

            /*Add offspring to population*/
            population.push_back(off1);
            costHeap->push_min(off1->getCost(), population.size() - 1);

            /*Add mutated to population*/
            population.push_back(off2);
            costHeap->push_min(off2->getCost(), population.size() - 1);

            /*If the size of population achieves $maxPopulation, then select survivors*/
            if(population.size() > maxPopulation) {
                printf("----- survivors selection\n");
                population = selectSurvivors(costHeap, population, sizePopulation, dataFrame);
                //bestSolution = population[costHeap->front_min().second];
            }

            /*If $itDiv iterations happened without improving the best solution, then diversify population*/ 
            if( ((it-lastImprovement) >= itDiv) && ((it-lastDiv) >= itDiv) ) {
                printf("----- diversification\n");
                lastDiv = it;
                population = diversifyPopulation(costHeap, population, sizePopulation, 2*sizePopulation, m, dataFrame, conflictGraph);
                if(costHeap->front_min().first < bestSolution->getCost()) {
                    lastImprovement = it;
                }
                bestSolution = population[costHeap->front_min().second];
            }

            it++;
            printf("%d) %.15f %s %d %s %d \n", it, bestSolution->getCost(), "size_pop =", (int)population.size(), "last_imp =", lastImprovement);
        }
        
        // Stop the running time clock
        double elapsedSecs = double(clock() - begin) / CLOCKS_PER_SEC;
        
        printSummary(bestSolution, elapsedSecs);

        /*double costVerified = verifyCost(bestSolution, dataFrame);
        printf("%.15f ", costVerified);*/
        
        delete [] annotated;
        delete dataFrame;
        for(int i = 0; i < population.size(); i++) {
            delete population[i];
        }
    }
}

/*Create a Solver instance*/
Solver* createSolver(DataFrame* dataFrame, int* solution) {
    Solver* solver;

    /*Create a Solver instance. The object to be created must be of a class that inherit Solver*/
    if(solverId == MEAN_ID) {
        solver = new KMeansSolver(dataFrame, solution, MEAN_ID);
    } else if(solverId == MEDIAN_ID) {
        solver = new KMediansSolver(dataFrame, solution, MEDIAN_ID);
    } else if(solverId == MEDOID_ID) {
        solver = new KMedoidsSolver(dataFrame, solution, MEDOID_ID);
    } else if(solverId == CSG_ID) {
        solver = new CsgSolver(dataFrame, solution, CSG_ID);
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
    
    //Call the main function, which will perform the genetic loop
    demo(1607, fileName, m);

    return 0;
}