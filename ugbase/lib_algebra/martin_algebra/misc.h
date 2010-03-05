/*
 *  misc.h
 *  flexamg
 *
 *  Created by Martin Rupp on 07.09.09.
 *  Copyright 2009 G-CSC. All rights reserved.
 *
 */
#pragma once
#include <iostream>
#include <iomanip>
#include <string>
#include <assert.h>
using namespace std;

namespace ug {
#define TRUE 1
#define FALSE 0



////////////////////////////////////////////////////////////////////////////////

#define boldredcolor "\x1b[1;31m"
#define boldgreencolor "\x1b[1;32m"
#define boldbluecolor "\x1b[1;34m"

#define redcolor "\x1b[0;31m"
#define greencolor "\x1b[0;32m"
#define bluecolor "\x1b[0;34m"

#define normalcolor "\x1b[0;0m"


////////////////////////////////////////////////////////////////////////////////
struct COMPILE_TIME_ERROR {};

template<bool> struct CompileTimeError;
template<> struct CompileTimeError<true> { };
#define CT_CHECK(expression, the_text)  \
	(CompileTimeError<(expression) != 0>() );


////////////////////////////////////////////////////////////////////////////////

//!
//! use this variable for 
static int never_happens = 0;
#define FORCE_CREATION if(never_happens)

extern int flexamg_dimensions;

#define OUTPUT_DIR "/Users/mrupp/matrices/"

#define IF_PRINTLEVEL(i) if(i <= 0)

////////////////////////////////////////////////////////////////////////////////
// PREFETCHING
//--------------
// 
//! gcc __builtin_prefetch(addr, readwrite =1 read = 0, locality)
//! locality=0..3, 0 dont leave in caches (default 3)
//! http://gcc.gnu.org/onlinedocs/gcc-4.1.2/gcc/Other-Builtins.html#Other-Builtins
#define prefetchRead(p) __builtin_prefetch(p, 0)
#define prefetchReadWrite(p) __builtin_prefetch(p, 1)
// note: p doesnt need to be a valid adress (for example, p can be NULL), but evaluating the expression (p) has to be possible

////////////////////////////////////////////////////////////////////////////////

string nrstring(double d);
string nrstring(int i);

////////////////////////////////////////////////////////////////////////////////
//#define DEBUG
#ifdef DEBUG

//#define IFDEBUG(a) a
//!
//! extended assert marco. one can insert a sequence passed to cout as second parameter
//! example: ASSERT2(nrOfElements < maxSize, "Number of Elements (" << nrOfElements+1 << ") exceeds arrays size (" << maxSize << ").");
//! void cArray::add(T &t) { ASSERT2(nrOfElements < maxSize, *this); , when ofstream operator <<  is supported by cArray.
//! ASSERT2 also prints out callstack information via backtrace.
#define ASSERT2(b, coutstuff) \
if(!(b)) \
{ \
cout << endl << redcolor << "===================== ERROR =====================" << normalcolor; \
cout << endl << endl << __func__ << ": " << coutstuff << endl;  \
cout << endl << "BACKTRACE:" << endl; \
print_trace(); \
cout << endl; \
cout.flush(); \
assert(b);  \
}

//!
//! print_trace: function to print callstack like backtrace in gdb.
#include <execinfo.h>
inline void print_trace ()
{
	void *array[256]; char **strings;	
	size_t size = backtrace (array, 256);
	strings = backtrace_symbols (array, size);	
	for (size_t i = 1; i < size; i++) // skip print_trace
		cout << strings[i] << endl;	
	free (strings);
}

#define ASSERT1(b) ASSERT2(b, "(no information available)");
#else

//#define IFDEBUG(a)
#define ASSERT2(a, b)
#define ASSERT1(a)
#endif

////////////////////////////////////////////////////////////////////////////////
inline double cut(double a, double e)
{
	if(a*a > e*e) return a;
	else return 0.0;
}



////////////////////////////////////////////////////////////////////////////////

//!
//! template struct for sorting some keys after values
//! for example, sorting a vector of ints and know original pos
template<typename T>
struct sortStruct
{
	int index; // for example "original" position.
	T sortValue;
	
	bool operator < (const sortStruct<T> &other) const
	{
		return sortValue < other.sortValue;
	}
};

enum Operation_type
{
	OPERATION_ADD, OPERATION_SET, OPERATION_SUB
};

inline Operation_type getAdd(Operation_type op)
{
	return ((op == OPERATION_SUB) ? OPERATION_SUB : OPERATION_ADD);
}

inline Operation_type getSub(Operation_type op)
{
	return ((op == OPERATION_SUB) ? OPERATION_ADD : OPERATION_SUB);
}

} // namespace ug

#include "positions.h"
#include "stopwatch.h"
