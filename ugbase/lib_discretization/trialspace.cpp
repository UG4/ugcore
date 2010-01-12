/*
 * trialspace.cpp
 *
 *  Created on: 12.05.2009
 *      Author: andreasvogel
 */

#include "trialspace.h"

namespace ug{

template<>
bool P1conform<Triangle>::evaluateShape(int nrShapeFct, const MathVector< RefDim >& locPos, number& value) const
{
	switch(nrShapeFct)
	{
		case 0: value = (1.-locPos.x-locPos.y); return true;
		case 1: value = locPos.x; return true;
		case 2: value = locPos.y; return true;
		default: return false;
	}
	return true;
};

template<>
bool P1conform<Triangle>::evaluateShape(int nrShapeFct, MathVector< RefDim > locPos[], number values[], int n) const
{
	for(int i = 0; i<n; i++)
	{
		switch(nrShapeFct)
		{
			case 0: values[i] = (1.-locPos[i].x-locPos[i].y); return true;
			case 1: values[i] = locPos[i].x; return true;
			case 2: values[i] = locPos[i].y; return true;
			default: return false;
		}
	}

	return true;
}

template<>
bool P1conform<Triangle>::evaluateShapeGrad(int nrShapeFct, const MathVector< RefDim >& locPos, ug::MathVector< RefDim >& value) const
{
	switch(nrShapeFct)
	{
	case 0: value[0] = -1.0;
			value[1] = -1.0; return true;
	case 1: value[0] = 	1.0;
			value[1] =  0.0; return true;
	case 2: value[0] =  0.0;
			value[1] =  1.0; return true;
	default: return false;
	}
	return true;
}

template<>
bool P1conform<Triangle>::positionOfDoF(int nrShapeFct, MathVector< RefDim >& value) const
{
	switch(nrShapeFct)
	{
	case 0: value[0] =  0.0;
			value[1] =  0.0; return true;
	case 1: value[0] = 	1.0;
			value[1] =  0.0; return true;
	case 2: value[0] =  0.0;
			value[1] =  1.0; return true;
	default: return false;
	}
	return true;
}


template<>
bool P1conform<Quadrilateral>::evaluateShape(int nrShapeFct, const MathVector< RefDim >& locPos, number& value) const
{
	switch(nrShapeFct)
	{
		 case 0: value = (1.-locPos.x)*(1.-locPos.y); return true;
		 case 1: value = locPos.x*(1.-locPos.y); return true;
		 case 2: value = locPos.x*locPos.y; return true;
		 case 3: value = (1.-locPos.x)*locPos.y; return true;
		 default: return false;
	}
	return true;
};

template<>
bool P1conform<Quadrilateral>::evaluateShape(int nrShapeFct, MathVector< RefDim > locPos[], number values[], int n) const
{
	for(int i = 0; i<n; i++)
	{
		switch(nrShapeFct)
		{
			case 0: values[i] = (1.-locPos[i].x)*(1.-locPos[i].y); return true;
			case 1: values[i] = locPos[i].x*(1.-locPos[i].y); return true;
			case 2: values[i] = locPos[i].x*locPos[i].y; return true;
			case 3: values[i] = (1.-locPos[i].x)*locPos[i].y; return true;
			default: return false;
		}
	}

	return true;
}

template<>
bool P1conform<Quadrilateral>::evaluateShapeGrad(int nrShapeFct, const MathVector< RefDim >& locPos, ug::MathVector< RefDim >& value) const
{
	switch(nrShapeFct)
	{
	case 0: value[0] = -1.0 + locPos.y;
			value[1] = -1.0 + locPos.x; return true;
	case 1: value[0] = 	1.0 - locPos.y;
			value[1] =      - locPos.x; return true;
	case 2: value[0] =        locPos.y;
			value[1] =        locPos.x; return true;
	case 3: value[0] =      - locPos.y;
			value[1] =  1.0 - locPos.x; return true;
	default: return false;
	}
	return true;
}

template<>
bool P1conform<Quadrilateral>::positionOfDoF(int nrShapeFct, MathVector< RefDim >& value) const
{
	switch(nrShapeFct)
	{
	case 0: value[0] =  0.0;
			value[1] =  0.0; return true;
	case 1: value[0] = 	1.0;
			value[1] =  0.0; return true;
	case 2: value[0] =  1.0;
			value[1] =  1.0; return true;
	case 3: value[0] =  0.0;
			value[1] =  1.0; return true;
	default: return false;
	}
	return true;
}


}
