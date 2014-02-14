/*
 * clock_fix.cpp
 *
 *  Created on: 14.02.2014
 *      Author: mrupp
 */




#include <time.h>
#include <unistd.h>
#include <stdio.h>

bool clock_fix_dummy()
{
	struct timespec tps;
	struct timespec res;
	return clock_gettime(CLOCK_REALTIME, &tps) && clock_getres(CLOCK_MONOTONIC, &res);

}
