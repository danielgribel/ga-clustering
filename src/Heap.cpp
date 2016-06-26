/************************************************************************************
Heap.cpp
Heap

Created by Daniel Gribel

This cpp file contains the Heap class definition.

A heap data structure to keep and manage the best solutions.
*************************************************************************************/
#include "Heap.h"

/*Heap constructor*/
HeapPdi::HeapPdi() {

}

/*Heap destructor*/
HeapPdi::~HeapPdi() {

}

/*Inserts an item considering a max heap*/
void HeapPdi::push_max(double cost, int val) {
    heap.push_back( std::pair<double, int>(cost, val) );
    push_heap(heap.begin(), heap.end());
}

/*Remove item with maximum value*/
int HeapPdi::pop_max() {
    make_heap(heap.begin(), heap.end());
    std::pair<double, int> pair = heap.front();
    pop_heap(heap.begin(), heap.end());
    heap.pop_back();
    return pair.second;
}

/*Inserts an item considering a min heap*/
void HeapPdi::push_min(double cost, int val) {
    heap.push_back( std::pair<double, int>(cost, val) );
    push_heap(heap.begin(), heap.end(), std::greater< std::pair<double, int> >());
}

/*Remove item with minimum value*/
int HeapPdi::pop_min() {
    make_heap(heap.begin(), heap.end(), std::greater< std::pair<double, int> >());
    std::pair<double, int> pair = heap.front();
    pop_heap(heap.begin(), heap.end(), std::greater< std::pair<double, int> >());
    heap.pop_back();
    return pair.second;
}

/*Gets the first item in the heap if we consider a max heap*/ 
std::pair<double, int> HeapPdi::front_max() {
    make_heap(heap.begin(), heap.end());
    return heap.front();
}

/*Gets the first item in the heap if we consider a min heap*/ 
std::pair<double, int> HeapPdi::front_min() {
    make_heap(heap.begin(), heap.end(), std::greater< std::pair<double, int> >());
    return heap.front();
}

/*Sets a new heap data*/
void HeapPdi::setHeap(std::vector< std::pair<double, int> > aHeap) {
    std::vector< std::pair<double, int> >(heap).swap(heap);
    heap = aHeap;
}

/*Clears the heap -- set it to be empty*/
void HeapPdi::clearHeap() {
    heap.resize(0);
    std::vector< std::pair<double, int> >(heap).swap(heap);
}