//*****************************************************************
//  LinkedList.h
//  HashTable
//
//  Created by Karlina Beringer on June 16, 2014.
//
//  This header file contains the Linked List class declaration.
//  Hash Table array elements consist of Linked List objects.
//*****************************************************************

#ifndef LinkedList_h
#define LinkedList_h

#include <iostream>
#include <string>
#include <algorithm>

using namespace std;

//*****************************************************************
// List items are keys with pointers to the next item.
//*****************************************************************
struct Item {
    double cost;
    int* cardinality;
    Item * next;
};

//*****************************************************************
// Linked lists store a variable number of items.
//*****************************************************************
class LinkedList
{
private:
    // Head is a reference to a list of data nodes.
    Item * head;
    
    // Length is the number of data nodes.
    int length;
    
public:
    // Constructs the empty linked list object.
    // Creates the head node and sets length to zero.
    LinkedList();
    
    // Inserts an item at the end of the list.
    void insertItem( Item * newItem );
    
    // Removes an item from the list by item key.
    // Returns true if the operation is successful.
    bool existItem(int* cardinality, double cost, int m);
    
    // Displays list contents to the console window.
    void printList();
    
    // Returns the length of the list.
    int getLength();
    
    // De-allocates list memory when the program terminates.
    ~LinkedList();
};

#endif

//*****************************************************************
// End of File
//*****************************************************************