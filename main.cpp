//
//  main.cpp
//  FiniteElementClass
//
//  Created by Victor Vidal Suriel on 3/18/16.
//  Copyright © 2016 Victor Vidal Suriel. All rights reserved.
//  test

#include <Accelerate/Accelerate.h>
#include <iostream>
#include <fstream>
#include "FEMclass.hpp"
using namespace std;

int main(int argc, const char * argv[]) {
/*    if (argc != 8) {
        cout << "ERROR: Put in 7 arguments\n";
        exit(-999);
    }
    double L=atof(argv[1]), a=atof(argv[2]), P=atof(argv[3]), Q=atof(argv[4]), f=atof(argv[5]), EI=atof(argv[6]); //Fixed parameters from problem defined by user
    int numEl = atoi(argv[7]);*/
    double L = 5, a = 2, P = -5, Q = -3, f = -1, EI = 1;
    int numEl = 60;    Element::setNumNodes(2);
    Element::setNumDOF(4);
    if (numEl % 3 != 0) {
        cout << "ERROR: number of elements must be divisible by 3\n";
        exit(-1);
    }
    int numDirDOF = 2;
    int dirDOF[numDirDOF];
    dirDOF[0] = 0;
    dirDOF[1] = 4;
    
    int numEssDOF = 2;
    int essDOF[numEssDOF];
    essDOF[0] = 2;
    essDOF[1] = 6;

    double essVal[numEssDOF];
    essVal[0] = P;
    essVal[1] = Q;

    int totNodes = numEl+1;
    double coords[totNodes];
    coords[0] = 0;
    coords[1] = L/2;
    coords[2] = L;
    coords[3] = L + a;
    
    int numSubEl = numEl/3;
    for (int i=1; i<numSubEl; i++) {
        coords[i+3] = i*L/(2*numSubEl);
        coords[i+numSubEl+2] = i*L/(2*numSubEl) + L/2;
        coords[i+2*numSubEl+1] = i*a/(numSubEl) + L;
    }
    
    Element * element = new Element[numEl];
    element[0].setNode(0, 0);
    element[0].setNode(1, 4);
    
    element[1].setNode(0, numSubEl + 2);
    element[1].setNode(1, 1);
    
    element[2].setNode(0, 1);
    element[2].setNode(1, numSubEl + 3);
    
    element[3].setNode(0, 2*numSubEl + 1);
    element[3].setNode(1, 2);
    
    element[4].setNode(0, 2);
    element[4].setNode(1, 2*numSubEl + 2);
    
    element[5].setNode(0, 3*numSubEl);
    element[5].setNode(1, 3);
    
    for (int i=1; i<numSubEl-1; i++) {
        element[i+5].setNode(0, i+3);
        element[i+5].setNode(1, i+4);
        
        element[i+numSubEl+3].setNode(0, i+numSubEl+2);
        element[i+numSubEl+3].setNode(1, i+numSubEl+3);
        
        element[i+2*numSubEl+1].setNode(0, i+2*numSubEl+1);
        element[i+2*numSubEl+1].setNode(1, i+2*numSubEl+2);
    }
    
    for (int i=0; i < numEl; i++) {
        element[i].setDOF(0, 2*element[i].getNode(0));
        element[i].setDOF(1, 2*element[i].getNode(0) + 1);
        element[i].setDOF(2, 2*element[i].getNode(1));
        element[i].setDOF(3, 2*element[i].getNode(1) + 1);
    }

    int totDOF = Element::getNumNodes()*totNodes;
    
    double * Stiff = new double[totDOF*totDOF];
    double * Load  = new double[totDOF];
    double elStiff[Element::getNumDOF()*Element::getNumDOF()];
    
    double h = 0;
    for (int i=0; i<numEl; i++) {
        h = coords[element[i].getNode(1)] - coords[element[i].getNode(0)];
        
        elStiff[0]  =  EI*12.0/pow(h, 3);
        elStiff[1]  =  EI*6.0/pow(h, 2);
        elStiff[2]  = -EI*12.0/pow(h, 3);
        elStiff[3]  =  EI*6.0/pow(h, 2);
        
        elStiff[4]  =  EI*6.0/pow(h, 2);
        elStiff[5]  =  EI*4.0/(h);
        elStiff[6]  = -EI*6.0/pow(h, 2);
        elStiff[7]  =  EI*2.0/(h);
    
        elStiff[8]  = -EI*12.0/pow(h, 3);
        elStiff[9]  = -EI*6.0/pow(h, 2);
        elStiff[10] =  EI*12.0/pow(h, 3);
        elStiff[11] = -EI*6.0/pow(h, 2);
        
        elStiff[12] =  EI*6.0/pow(h, 2);
        elStiff[13] =  EI*2.0/(h);
        elStiff[14] = -EI*6.0/pow(h, 2);
        elStiff[15] =  EI*4.0/(h);
        
        element[i].updateStiff(Stiff, elStiff, totDOF);
        
        Load[element[i].getDOF(0)] +=  f*(h)/2.0;
        Load[element[i].getDOF(1)] +=  f*pow(h,2)/12.0;
        Load[element[i].getDOF(2)] +=  f*(h)/2.0;
        Load[element[i].getDOF(3)] += -f*pow(h,2)/12.0;
    }
    double rawStiff[numDirDOF][totDOF];
    double rawLoad[numDirDOF];
    
    for (int i=0; i<numDirDOF; i++) {
        for (int j=0; j<totDOF; j++) {
            rawStiff[i][j] = Stiff[totDOF*j + dirDOF[i]];
            Stiff[totDOF*j + dirDOF[i]] = 0;
            Stiff[totDOF*dirDOF[i] + j] = 1;
        }
        Stiff[dirDOF[i]*totDOF + dirDOF[i]] = 1;
        rawLoad[i] = Load[dirDOF[i]];
        Load[dirDOF[i]] = 0;
    }
    for (int i=0; i<numEssDOF; i++)
        Load[essDOF[i]] += essVal[i];
    
    double * x = new double[totDOF];// This stores the value of the solution
    for(int i=0; i<totDOF; i++) {
        x[i]=Load[i];
    }
    
    __CLPK_integer n, lda, ldb, nrhs, info;
    n = lda = ldb = totDOF;
    nrhs=1;
    __CLPK_integer * ipiv = new __CLPK_integer[totDOF];
    dgesv_(&n, &nrhs, Stiff, &lda, ipiv, x, &ldb, &info);

    double dirForce[numDirDOF];//Stores the residual force at points A and C
    
    for (int i=0; i<numDirDOF; i++) {
        for (int j=0; j<totDOF; j++) {
            dirForce[i] += rawStiff[i][j]*x[j];
        }
        dirForce[i] -= rawLoad[i];
    }
    //Print the forces at A and C onto the screen
    for (int i=0; i<numDirDOF; i++)
        cout << "The force at point " << i << " is " << dirForce[i] << endl;
    
    //Writes the values of coords and the even values of x (starting at zero), which have the values of the displacement. These are written onto a text file for plotting by MATLAB
    ofstream outFile;
    outFile.open("FEM2.txt");
    for (int i=0; i<totNodes; i++) {
        outFile << coords[i] << " " << x[2*i] << endl;
    }
    outFile.close();
        return 0;
}
