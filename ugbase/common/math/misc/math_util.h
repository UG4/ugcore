//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d14

#ifndef __H__UGMATH__MATH_UTIL__
#define __H__UGMATH__MATH_UTIL__

#include <cstdlib>

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	urand
///	uniform distributed random numbers in [lowerBound, upperBound[. Use srand to set a seed.
template <class TNumber>
TNumber
urand(TNumber lowerBound, TNumber upperBound)
{
	long t = rand();
	if(t == RAND_MAX)
		t -= 1;

	return lowerBound + (TNumber)((upperBound - lowerBound) * ((float)t / (float)RAND_MAX));
}

////////////////////////////////////////////////////////////////////////
//	clip
///	clips a number to the given interval [lowerBound, upperBound].
template <class TNumber>
TNumber
clip(TNumber val, TNumber lowerBound, TNumber upperBound)
{
	if(val > upperBound)
		return upperBound;
	else if(val < lowerBound)
		return lowerBound;
	return val;
}

}//	end of namespace

#endif
