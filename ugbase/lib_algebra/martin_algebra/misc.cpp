/*
 *  misc.cpp
 *  flexamg
 *
 *  Created by Martin Rupp on 07.09.09.
 *  Copyright 2009 G-CSC. All rights reserved.
 *
 */



#include "misc.h"

const char *boldredcolor = "\x1b[1;31m";
const char *boldgreencolor = "\x1b[1;32m";
const char *boldbluecolor = "\x1b[1;34m";

const char *redcolor = "\x1b[0;31m";
const char *greencolor = "\x1b[0;32m";
const char *bluecolor = "\x1b[0;34m";

const char *normalcolor = "\x1b[0;0m";

string nrstring(double d)
{
	char s[255];
	sprintf(s, "%g", d);
	return string(s);
}

string nrstring(int i)
{
	char s[255];
	sprintf(s, "%d", i);
	return string(s);
}

void spaceout(int n)
{
	for(int i=0;i<n; i++)
		cout << " ";
}

int *parentIndex[32];
#ifdef DEBUG
void print_trace ()
{
	void *array[256];
	size_t size;
	char **strings;
	size_t i;
	
	size = backtrace (array, 256);
	strings = backtrace_symbols (array, size);
	
	for (i = 1; i < size; i++)
		cout << strings[i] << endl;
	
	free (strings);
}

void my_break_func()
{
	print_trace();
}
#endif