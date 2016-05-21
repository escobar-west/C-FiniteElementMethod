//
//  FEMclass.hpp
//  FiniteElementClass
//
//  Created by Victor Vidal Suriel on 3/18/16.
//  Copyright © 2016 Victor Vidal Suriel. All rights reserved.
//

#ifndef FEMclass_hpp
#define FEMclass_hpp

class Element {
private:
    static int numNodes; //Number of nodes each element has
    static int numDOF;   //Number of DOF each element has
    int * Node; //points to array of Nodes
    int * DOF;  //points to array of DOF
public:
    Element();
    Element(int * new_nodes, int * new_DOF);
    Element(const Element &object);
    Element& operator=(const Element &object);
    ~Element();
    static void setNumNodes(int new_number);
    static void setNumDOF(int new_number);
    void setNode(int pos, int val);
    void setDOF(int pos, int val);
    static int getNumNodes() {return numNodes;}
    static int getNumDOF()   {return numDOF;}
    int getNode(int pos);
    int getDOF(int pos);
    void updateStiff(double * A, double * a, int totDOF);
};
#endif /* FEMclass_hpp */
