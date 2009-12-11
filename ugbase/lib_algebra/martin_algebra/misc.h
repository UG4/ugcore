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

#define INTBITSIZE (sizeof(int)*8)

extern const char *boldredcolor, *boldgreencolor, *boldbluecolor, *redcolor, *greencolor, *bluecolor, *normalcolor ;

#define IF_PRINTLEVEL(i) if(i <= 0)

// gcc __builtin_prefetch(addr, readwrite =1 read = 0, locality)
// locality=0..3, 0 dont leave in caches (default 3)
// http://gcc.gnu.org/onlinedocs/gcc-4.1.2/gcc/Other-Builtins.html#Other-Builtins
#define prefetchRead(p) __builtin_prefetch(p, 0)
#define prefetchReadWrite(p) __builtin_prefetch(p, 1)

void spaceout(int n);
#ifdef DEBUG
#include <execinfo.h>
void print_trace();
void my_break_func();
#define BREAKON(a) if(a) my_break_func();
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
#define ASSERT(b) ASSERT2(b, "(no information available)");
#else
#define BREAKON(a)
#define ASSERT2(a, b)
#define ASSERT(a)
#endif



#define UNROLL 32
#ifdef UNROLL
#define FOR_UNROLL_FWD(i, for_start, for_length, unroll_length, _ex_) \
{ \
int i; \
for(i=for_start; i < for_length%unroll_length; i++) {_ex_; } \
while(i<for_length)\
	for(int for_unroll_temp2=0; for_unroll_temp2 < unroll_length; for_unroll_temp2++, i++) {_ex_;}  \
}
#else
#define FOR_UNROLL_FWD(i, for_start, for_length, unroll_length, _ex_) \
{ \
for(int i=for_start; i < for_length; i++) { _ex_; } \
}
#endif


inline double cut(double a, double e)
{
	if(a*a > e*e) return a;
	else return 0.0;
}

// oder Vector<bool>
class bitmarker
{
public:
	inline bitmarker(int _length)
	{
		length = _length;
		bitarray = new int[length / INTBITSIZE +1];
		memset(bitarray, 0, sizeof(int)*(length/INTBITSIZE +1));
	}
	
	inline ~bitmarker()
	{
		delete [] bitarray;
	}
	
	inline bool getbit(int i) const
	{
		return (bool) (bitarray[i/INTBITSIZE] & (1 << (i % INTBITSIZE)));
	}
		
	inline void setbit(int i, int b)
	{
		if(b)
			bitarray[i/INTBITSIZE] |= (1 << (i % INTBITSIZE));
		else
		{
			int x = ~((int)(1 << (i % INTBITSIZE)));
			bitarray[i/INTBITSIZE] &= x;
		}
	}
	
	void print() const
	{
		for(int i=0; i<length; i++)
			cout << getbit(i) ? "1" : "0";
		cout << endl;
	}
	
private:
	int length;
	int *bitarray;
	
};

struct pos2d
{
	double x, y;
};

pos2d GetPosForIndex(int i);

/*
class cBigmem
{
public:
	cBigmem(int bytes)
	{
		bigmem = new char[bytes];
	}
	
	~cBigmem()
	{
		delete[] bigmem;
	}
	
	bool lock()
	{
		ASSERT(!memlocked);
		memlocked = 1;
	}
	bool unlock()
	{
		ASSERT(memlocked);
		memlocked = 0;
	}

	void *getmem(int bytes)
	{
		ASSERT(endmem - freemem >= bytes);
		void *a = (void*) freemem;
		freemem += bytes;
		return a;
	}
		
	void *getmem(int nrofitems, int sizeofitem)
	{
		return getmem(nrofitems*sizeofitem);
	}	
	
private:
	char *bigmem;
	char *freemem;
	char *endmem;
	bool memlocked;	
};

extern cBigmem bigmem;
*/


void writePosToStream(ostream &out);

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

class stopwatch
{
public:
	stopwatch() 
	{
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
