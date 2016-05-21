//
//  FEMclass.cpp
//  FiniteElementClass
//
//  Created by Victor Vidal Suriel on 3/18/16.
//  Copyright © 2016 Victor Vidal Suriel. All rights reserved.
//

#include <iostream>
#include "FEMclass.hpp"
using namespace std;

int Element::numNodes = -1;
int Element::numDOF   = -1;

Element::Element() {
    Node = new int[numNodes];
    DOF = new int[numDOF];
}
    
Element::Element(int * new_nodes, int * new_DOF) {
    Node = new int[numNodes];
    for (int i=0; i < numNodes; i++)
        Node[i] = new_nodes[i];
    
    DOF = new int[numDOF];
    for (int i=0; i < numDOF; i++)
        DOF[i] = new_DOF[i];
}

Element::Element(const Element &object) {
    Node = new int[numNodes];
    for (int i=0; i < numNodes; i++)
        Node[i] = object.Node[i];
    
    DOF = new int[numDOF];
    for (int i=0; i < numDOF; i++)
        DOF[i] = object.DOF[i];
}

Element& Element::operator = (const Element & object) {
    if (this != &object) {
        if (Node != NULL)
            delete[] Node;
        
        if (DOF != NULL)
        delete[] DOF;
        
        Node = new int[numNodes];
        for (int i=0; i < numNodes; i++)
            Node[i] = object.Node[i];
        
        DOF = new int[numDOF];
        for (int i=0; i < numDOF; i++)
            DOF[i] = object.DOF[i];
    }
    return *this;
}

Element::~Element() {
    delete[] Node;
    delete[] DOF;
}

void Element::setNumNodes(int new_number) {
    if( numNodes != -1) {
        cout << "ERROR: attempted to change the number of element nodes twice\n";
        exit(-1);
    }
    if( new_number<=0 ) {
        cout << "ERROR: number of nodes must be positive\n";
        exit(-1);
    }
    numNodes = new_number;
}

void Element::setNumDOF(int new_number) {
if( numDOF != -1) {
    cout << "ERROR: attempted to change the number of element DOF twice\n";
    exit(-1);
    }
    if( new_number<0 ) {
        cout << "ERROR: cannot set number of element DOF to be negative\n";
        exit(-1);
    }
    numDOF = new_number;
    }

void Element::setNode(int pos, int val) {
    if (pos < 0 || numNodes <= pos) {
        cout << "ERROR: attempted to set non-existant node\n";
        exit(-1);
    }
    Node[pos] = val;
}

void Element::setDOF(int pos, int val) {
    if (pos < 0 || numDOF <= pos) {
        cout << "ERROR: attempted to set non-existant DOF\n";
        exit(-1);
    }
    DOF[pos] = val;
}

int Element::getNode(int pos) {
    if (pos < 0 || numNodes <= pos) {
        cout << "ERROR: attempted to get non-existant node\n";
        exit(-1);
    }
    return Node[pos];
}

int Element::getDOF(int pos) {
    if (pos < 0 || numDOF <= pos) {
        cout << "ERROR: attempted to get non-existant DOF\n";
        exit(-1);
    }
    return DOF[pos];
}

void Element::updateStiff(double * A, double * a, int totDOF) {
    for (int i=0; i< numDOF; i++)
        for (int j=0; j<numDOF; j++)
            A[totDOF*getDOF(j) + getDOF(i)] += a[numDOF*j + i];
}

/*Just stick this code in wherever you need to invert a matrix

void BLASsolver(double * A, double * b, int totDOF) {
__CLPK_integer n, lda, ldb, nrhs, info;
n = lda = ldb = totDOF;
nrhs=1;
__CLPK_integer * ipiv = new __CLPK_integer[totDOF];
dgesv_(&n, &nrhs, A, &lda, ipiv, b, &ldb, &info);
}
*/
