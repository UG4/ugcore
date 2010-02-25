/*
 *  algebra.h
 *  flexamg
 *
 *  Created by Martin Rupp on 05.01.10.
 *  Copyright 2010 . All rights reserved.
 *
 */

#include <iostream>
using namespace std;
#include <math.h>

#include "misc.h"

//#include "SparseMatrix.h"

#include "amg.h"
#include "famg.h"
#include "linearsolver.h"



#define EASY_MATRIX // set h = 1.0
#define NINE_POINT

#ifndef NINE_POINT
#define AGGRESSIVE_COARSENING
#endif
#define UNKNOWN_NR 1
//#define ONE_LEVEL

#include "arrayStorage.h"

#if UNKNOWN_NR > 1

//typedef fixedStorage myStorageType;
typedef variableStorage myStorageType;

typedef blockDenseMatrix<double, myStorageType, UNKNOWN_NR, UNKNOWN_NR > myBlockMat;

typedef blockVector<double, myStorageType, UNKNOWN_NR> myBlockVec;

#else

typedef double myBlockMat;
typedef double myBlockVec;

#endif


typedef SparseMatrix<myBlockMat> myMatrix;
typedef Vector<myBlockVec> myVector;