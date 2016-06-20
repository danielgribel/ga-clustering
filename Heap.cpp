#include "Heap.h"

HeapPdi::HeapPdi() {

}

HeapPdi::~HeapPdi() {

}

void HeapPdi::push_max(double cost, int val) {
    heap.push_back( std::pair<double, int>(cost, val) );
    push_heap(heap.begin(), heap.end());
}

int HeapPdi::pop_max() {
    make_heap(heap.begin(), heap.end());

    std::pair<double, int> pair = heap.front();

    //int val = heap.front();
    pop_heap(heap.begin(), heap.end());
    heap.pop_back();
    return pair.second;
}

void HeapPdi::push_min(double cost, int val) {
    heap.push_back( std::pair<double, int>(cost, val) );
    push_heap(heap.begin(), heap.end(), std::greater< std::pair<double, int> >());
}

int HeapPdi::pop_min() {
    make_heap(heap.begin(), heap.end(), std::greater< std::pair<double, int> >());
    
    std::pair<double, int> pair = heap.front();
    //int val = heap.front();
    pop_heap(heap.begin(), heap.end(), std::greater< std::pair<double, int> >());
    heap.pop_back();
    return pair.second;
}

std::pair<double, int> HeapPdi::front_max() {
    make_heap(heap.begin(), heap.end());
    //int val = heap.front();
    return heap.front();
}

std::pair<double, int> HeapPdi::front_min() {
    make_heap(heap.begin(), heap.end(), std::greater< std::pair<double, int> >());
    //int val = heap.front();
    return heap.front();
}

void HeapPdi::setHeap(std::vector< std::pair<double, int> > aHeap) {
    std::vector< std::pair<double, int> >(heap).swap(heap);
    heap = aHeap;
}

void HeapPdi::clearHeap() {
    heap.resize(0);
    std::vector< std::pair<double, int> >(heap).swap(heap);
}