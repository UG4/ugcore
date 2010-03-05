/*
 *  algebra.h
 *  flexamg
 *
 *  Created by Martin Rupp on 05.01.10.
 *  Copyright 2010 G-CSC, University of Frankfurt. All rights reserved.
 *
 */
#pragma once

#include <iostream>
using namespace std;

#include <math.h>

#include "misc.h"

//#include "SparseMatrix.h"

#include "amg.h"
//#include "famg.h"
#include "linearsolver.h"



//#define EASY_MATRIX // set h = 1.0

#ifndef UNKNOWN_NR
#define UNKNOWN_NR 1
#endif

//#define VARIABLE_STORAGE

namespace ug{
#if UNKNOWN_NR > 1
	
	typedef SparseMatrix<double> myMatrix2;
	
#ifdef VARIABLE_STORAGE
	typedef variableStorage myStorageType;
#else
	typedef fixedStorage myStorageType;
#endif
	
	typedef blockDenseMatrix<double, myStorageType, UNKNOWN_NR, UNKNOWN_NR > myBlockMat;
	
	typedef blockVector<double, myStorageType, UNKNOWN_NR> myBlockVec;
		
#else
	
	typedef double myBlockMat;
	typedef double myBlockVec;	
#endif	
	
	typedef SparseMatrix<myBlockMat> myMatrix;
	typedef Vector<myBlockVec> myVector;
}