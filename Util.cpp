#include "Util.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <queue>
#include <math.h>

using namespace std;

void shuffle(int *myArray, size_t n) {
    for(int i = 0; i < n; i++) {
        myArray[i] = i;
    }
    if(n > 1) {
        size_t i;
        for (int i = 0; i < n - 1; i++) {
            size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
            int t = myArray[j];
            myArray[j] = myArray[i];
            myArray[i] = t;
        }
    }
}

int* getRandomSolution(int n, int m) {
    int* s = new int[n];
    
    for(int i = 0; i < n; i++) {
        s[i] = rand() % m + 0;
    }

    return s;
}

int* perturbation(int* s0, int k, int n, int m) {
    int* s = new int[n];
    int randItems[n];

    for(int i = 0; i < n; i++) {
        s[i] = s0[i];
    }

    shuffle(randItems, n);

    for(int i = 0; i < k; i++) {
        s[randItems[i]] = rand() % m + 0;
    }

    return s;
}

double sumSquare(double* a, double* b, int d) {
    double sum = 0.0;
    double diff;
    for(int i = 0; i < d; i++) {
        diff = a[i] - b[i];
        sum = sum + diff*diff;
    }
    return sum;
}

double euclidean(double* a, double* b, int d) {
    double sum = sqrt(sumSquare(a, b, d));
    return sum;
}

double manhattan(double* a, double* b, int d) {
    double sum = 0.0;
    for(int i = 0; i < d; i++) {
        sum = sum + abs(a[i] - b[i]);
    }
    return sum;
}

double getDistance(double* a, double* b, int d) {
    return euclidean(a, b, d); 
}

double* kLargest(double* myArray, int m, int k) {
    double* A = new double[m];
    for(int i = 0; i < m; i++) {
        A[i] = myArray[i];
    }

    double* largests = new double[k];

    nth_element(A, A + (m-k), A + m);
    int j = 0;

    for(int i = m-k; i < m; i++) {
        largests[j] = A[i];
        j++;
    }

    delete [] A;
    
    return largests;
}

int* randAnnotated(int n, int m) {
    int in, im;
    int* annotated = new int[m];
    int rn, rm;
    im = 0;

    for(in = 0; in < n && im < m; ++in) {
        rn = n - in;
        rm = m - im;
        if(rand() % rn < rm) {
            annotated[im++] = in;
        }
    }
    //assert(im == m);

    return annotated;
}

vector<int> kSmallestIndices(vector<double> v, int k) {
    priority_queue<pair< double, int>, vector < pair< double, int> >, smallerOp > q;
    
    for(int i = 0; i < v.size(); ++i) {
        q.push(pair<double, int>(v[i], i));
    }
    
    vector<int> n_closest(k);
    int ki;
    
    for(int i = 0; i < k; ++i) {
        ki = q.top().second;
        n_closest[i] = ki;
        q.pop();
    }

    return n_closest;
}

vector<int> kSmallestIndices(double* v, int n, int k) {
    priority_queue<pair< double, int>, vector < pair< double, int> >, smallerOp > q;

    for(int i = 0; i < n; ++i) {
        q.push(pair<double, int>(v[i], i));
    }
    
    vector<int> n_closest(k);
    int ki;
    
    for(int i = 0; i < k; ++i) {
        ki = q.top().second;
        n_closest[i] = ki;
        q.pop();
    }

    return n_closest;
}

/*int main() {
    vector<double> v = {5.2, 1.0, 0.01, 3.0, 0.002, -1.0, 20};
    int k = 3; // number of indices we need

    vector<int> n_closest = kSmallestIndices(v, k);

    for(int i = 0; i < k; ++i) {
        cout << n_closest[i] << endl;
    }

    return 0;
}*/