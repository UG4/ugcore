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
#include <assert.h>
using namespace std;
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

static int never_happens = 0;

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

void spaceout(int n);
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
inline double dabs(double a) { return a > 0 ? a : -a; }



////////////////////////////////////////////////////////////////////////////////
// stuff for position output
struct postype
{
	double x, y, z;
	friend std::ostream &operator << (std::ostream &out, const postype &p)
	{
		if(flexamg_dimensions == 2)
			out << p.x << " " << p.y;
		else
			out << p.x << " " << p.y << " " << p.z;
		return out;
	}
};


postype GetPosForIndex(int i);

void writePosToStream(ostream &out);
void writeToPosFile(const char *filename);

extern int *parentIndex[32];
static int GetOriginalIndex(int level, int i)
{
	while(level > 0)
		i = parentIndex[level--][i];
	return i;
}

static postype GetPosForIndexAtLevel(int i, int level)
{
	return GetPosForIndex(GetOriginalIndex(level, i));
}

////////////////////////////////////////////////////////////////////////////////
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
	friend ostream &operator << (ostream &out, stopwatch &s)
	{
		out << s.getTimeDiffMS() << " ms";
		return out;
	}
	
	double getTimeDiffMS()
	{
		if(bRunning) end = clock();
		return (clock()-beg)/((double)0.001*CLOCKS_PER_SEC);
	}
	
	void printTimeDiff()
	{
		cout << "took " << getTimeDiffMS() << " ms" << endl;
		cout.flush();
	}
private:
	clock_t beg, end;
	bool bRunning;
};


