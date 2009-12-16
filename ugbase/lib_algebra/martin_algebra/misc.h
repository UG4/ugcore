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
#include <string>

using namespace std;
#define TRUE 1
#define FALSE 0

string nrstring(double d);
string nrstring(int i);

extern const char *boldredcolor, *boldgreencolor, *boldbluecolor, *redcolor, *greencolor, *bluecolor, *normalcolor ;

#define IF_PRINTLEVEL(i) if(i <= 0)

// PREFETCHING
//--------------
// 
// gcc __builtin_prefetch(addr, readwrite =1 read = 0, locality)
// locality=0..3, 0 dont leave in caches (default 3)
// http://gcc.gnu.org/onlinedocs/gcc-4.1.2/gcc/Other-Builtins.html#Other-Builtins
#define prefetchRead(p) __builtin_prefetch(p, 0)
#define prefetchReadWrite(p) __builtin_prefetch(p, 1)
// note: p doesnt need to be a valid adress (for example, p can be NULL), but evaluating the expression (p) has to be possible

void spaceout(int n);

#ifdef DEBUG


//!
//! extended assert marco. one can insert a sequence passed to cout as second parameter
//! example: ASSERT2(nrOfElements < maxSize, "Number of Elements (" << nrOfElements+1 << ") exceeds arrays size (" << maxSize << ").");
//! void cArray::add(T &t) { ASSERT2(nrOfElements < maxSize, *this); , when ofstream operator <<  is supported by cArray.
//! ASSERT2 also prints out callstack information via backtrace.
#define ASSERT2(b, coutstuff) \
if(!(b)) \
{ \
cout << endl << redcolor << "==================================================== ERROR ====================================================" << normalcolor; \
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

#define ASSERT(b) ASSERT2(b, "(no information available)");
#else

#define ASSERT2(a, b)
#define ASSERT(a)
#endif


inline double cut(double a, double e)
{
	if(a*a > e*e) return a;
	else return 0.0;
}



struct pos2d
{
	double x, y;
};
extern pos2d *positions;
extern int iNrOfPositions;

pos2d GetPosForIndex(int i);

void writePosToStream(ostream &out);
void writeToPosFile(const char *filename);

extern int *parentIndex[32];
static int GetOriginalIndex(int level, int i)
{
	while(level > 0)
		i = parentIndex[level--][i];
	return i;
}

static pos2d GetPosForIndexAtLevel(int i, int level)
{
	return GetPosForIndex(GetOriginalIndex(level, i));
}

//!
//! stopwatch class for quickly taking times
//! seems to be ok for measuring times > 100 ms
class stopwatch
{
public:
	stopwatch() 
	{
		// you cant be really sure when constructor is called
		beg = clock();
		bRunning = false;
	}
	void start()
	{
		cout.flush();
		beg = clock();
		bRunning = true;
	}
	void stop()
	{
		end = clock();
		bRunning = false;
	}
	
	double gettimediffMS()
	{
		if(bRunning) end = clock();
		return (clock()-beg)/((double)0.001*CLOCKS_PER_SEC);
	}
	
	void printTimeDiff()
	{
		cout << "took " << gettimediffMS() << " ms" << endl;
		cout.flush();
	}
private:
	clock_t beg, end;
	bool bRunning;
};

inline double dabs(double a) { return a > 0 ? a : -a; }
