//*****************************************************************
//  LinkedList.cpp
//  HashTable
//
//  Created by Karlina Beringer on June 16, 2014.
//
//  This header file contains the Linked List class declaration.
//  Hash Table array elements consist of Linked List objects.
//*****************************************************************

#include "LinkedList.h"

// Constructs the empty linked list object.
// Creates the head node and sets length to zero.
LinkedList::LinkedList() {
    head = new Item;
    head -> next = NULL;
    length = 0;
}

// Inserts an item at the end of the list.
void LinkedList::insertItem( Item * newItem ) {
    if (!head -> next)
    {
        head -> next = newItem;
        length++;
        return;
    }
    Item * p = head;
    Item * q = head;
    while (q)
    {
        p = q;
        q = p -> next;
    }
    p -> next = newItem;
    newItem -> next = NULL;
    length++;
}

bool checkSize(Item * p, int* cardinality, int m) {
    for(int i = 0 ; i < m; i++) {
        if(p->cardinality[i] != cardinality[i]) {
            return false;
        }
    }

    return true;
}

bool LinkedList::existItem(int* cardinality, double cost, int m) {
    Item * p = head;
    Item * q = head;
    
    while(q) {
        p = q;
        if ((p != head) && (checkSize(p, cardinality, m)) && (cost > p->cost - 0.00000001) && (cost < p->cost + 0.00000001)) {
            return true;
        }
        q = p->next;
    }
    return false;
}

// Displays list contents to the console window.
void LinkedList::printList()
{
    if (length == 0)
    {
        cout << "\n{ }\n";
        return;
    }
    Item * p = head;
    Item * q = head;
    cout << "\n{ ";
    while (q)
    {
        p = q;
        if (p != head)
        {
            if (p -> next) cout << ", ";
            else cout << " ";
        }
        q = p -> next;
    }
    cout << "}\n";
}

// Returns the length of the list.
int LinkedList::getLength() {
    return length;
}

// De-allocates list memory when the program terminates.
LinkedList::~LinkedList() {
    Item * p = head;
    Item * q = head;
    while (q)
    {
        p = q;
        q = p -> next;
        if (q) delete p;
    }
}