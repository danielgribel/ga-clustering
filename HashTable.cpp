//*****************************************************************
//  HashTable.cpp
//  HashTable
//
//  Created by Kar Beringer on June 18, 2014.
//
//  This header file contains the Hash Table class definition.
//  Hash Table array elements consist of Linked List objects.
//*****************************************************************

#include "HashTable.h"

// Constructs the empty Hash Table object.
HashTable::HashTable( int tableLength ) {
    if (tableLength <= 0) {
        tableLength = 211;
    } 
    array = new LinkedList[tableLength];
    length = tableLength;
}

// Returns an array location for a given item key.
int HashTable::hash(int* cardinality, int m) {
    sort(cardinality, cardinality + m);
    int value = 0;

    for(int i = 0; i < m; i++) {
        value = value + (i+1)*cardinality[i];
    }

    return value % length;
}

// Adds an item to the Hash Table.
void HashTable::insertItem( Item * newItem, int m ) {
    int index = hash(newItem -> cardinality, m);
    array[ index ].insertItem( newItem );
}

bool HashTable::existItem(int* cardinality, double cost, int m) {
    int index = hash(cardinality, m);
    return array[index].existItem(cardinality, cost, m);
}

// Display the contents of the Hash Table to console window.
void HashTable::printTable() {
    cout << "\n\nHash Table:\n";
    for ( int i = 0; i < length; i++ ) {
        cout << "Bucket " << i + 1 << ": ";
        array[i].printList();
    }
}

// Prints a histogram illustrating the Item distribution.
void HashTable::printHistogram() {
    cout << "\n\nHash Table Contains ";
    cout << getNumberOfItems() << " Items total\n";
    for ( int i = 0; i < length; i++ )
    {
        cout << i + 1 << ":\t";
        for ( int j = 0; j < array[i].getLength(); j++ )
            cout << " X";
        cout << "\n";
    }
}

// Returns the number of locations in the Hash Table.
int HashTable::getLength() {
    return length;
}

// Returns the number of Items in the Hash Table.
int HashTable::getNumberOfItems() {
    int itemCount = 0;
    for ( int i = 0; i < length; i++ )
    {
        itemCount += array[i].getLength();
    }
    return itemCount;
}

// De-allocates all memory used for the Hash Table.
HashTable::~HashTable() {
    delete [] array;
}