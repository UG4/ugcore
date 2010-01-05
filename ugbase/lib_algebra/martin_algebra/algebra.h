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

const double sigma = 0.3;
const double theta = 0.3;

#include "amg.h"
#include "ls.h"



#define EASY_MATRIX // set h = 1.0
//#define NINE_POINT

#ifndef NINE_POINT
#define AGGRESSIVE_COARSENING
#endif
#define UNKNOWN_NR 2
//#define ONE_LEVEL

#include "arrayStorage.h"

#if UNKNOWN_NR > 1
typedef fixedStorage myStorageType;
//typedef variableStorage myStorageType;
typedef fixedMatrix<myStorageType, UNKNOWN_NR, UNKNOWN_NR > MAT_TYPE;
typedef fixedVector<myStorageType, UNKNOWN_NR> VEC_TYPE;
#else
typedef double VEC_TYPE;
typedef double MAT_TYPE;

//typedef fixedMatrix<1> MAT_TYPE;
//typedef fixedVector<1> VEC_TYPE;

#endif

