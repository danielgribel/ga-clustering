/************************************************************************************
Heap.cpp
Heap

Created by Daniel Gribel

This header file contains the Heap class declaration.
*************************************************************************************/

#ifndef Heap_Pdi_H
#define Heap_Pdi_H

#include <algorithm>
#include <iostream>
#include <vector>
#include <iterator>
#include <functional>

class HeapPdi {
    private:
        /*Heap data structure. Each item in the heap is a pair <double, int>
        The double value indicates the cost of a solution
        The int value is the id of a solution*/
        std::vector< std::pair<double, int> > heap;

    public:
        /*Heap constructor*/
        HeapPdi();

        /*Heap destructor*/
        ~HeapPdi();
        
        /*Inserts an item considering a max heap*/
        void push_max(double cost, int val);
        
        /*Remove item with maximum value*/
        int pop_max();
        
        /*Inserts an item considering a min heap*/
        void push_min(double cost, int val);
        
        /*Remove item with minimum value*/
        int pop_min();
        
        /*Returns the heap data structure*/
        std::vector< std::pair<double, int> > getHeap() { return heap; };

        /*Gets the first item in the heap if we consider a max heap*/ 
        std::pair<double, int> front_max();
        
        /*Gets the first item in the heap if we consider a min heap*/
        std::pair<double, int> front_min();
                
        /*Sets a new heap data*/
        void setHeap(std::vector< std::pair<double, int> > aHeap);
        
        /*Clears the heap -- set it to be empty*/
        void clearHeap();
};

#endif