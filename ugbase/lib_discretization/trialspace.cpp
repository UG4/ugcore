/*
 * trialspace.cpp
 *
 *  Created on: 12.05.2009
 *      Author: andreasvogel
 */

#include "trialspace.h"

namespace ug{

bool operator==(const ElementDoFPattern& lhs, const ElementDoFPattern& rhs)
{
	for(uint i = 0; i < NUM_GEOMETRIC_BASE_OBJECTS; ++i)
	{
		if(lhs.m_dof_pattern[i] != rhs.m_dof_pattern[i]) return false;
	}
	return true;
}

// definition of static variables
std::map<TrialSpaceType, bool> TrialSpaces::m_bIsAdaptive = std::map<TrialSpaceType, bool>();
std::map<TrialSpaceType, std::map<int, ElementDoFPattern> > TrialSpaces::m_dimDoFPattern = std::map<TrialSpaceType, std::map<int, ElementDoFPattern> >();


template<>
bool P1conform<Triangle>::evaluate_shape(int nrShapeFct, const MathVector< RefDim >& locPos, number& value) const
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
bool P1conform<Triangle>::evaluate_shape(int nrShapeFct, MathVector< RefDim > locPos[], number values[], int n) const
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
bool P1conform<Triangle>::evaluate_shape_grad(int nrShapeFct, const MathVector< RefDim >& locPos, ug::MathVector< RefDim >& value) const
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
bool P1conform<Triangle>::position_of_dof(int nrShapeFct, MathVector< RefDim >& value) const
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
bool P1conform<Quadrilateral>::evaluate_shape(int nrShapeFct, const MathVector< RefDim >& locPos, number& value) const
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
bool P1conform<Quadrilateral>::evaluate_shape(int nrShapeFct, MathVector< RefDim > locPos[], number values[], int n) const
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
bool P1conform<Quadrilateral>::evaluate_shape_grad(int nrShapeFct, const MathVector< RefDim >& locPos, ug::MathVector< RefDim >& value) const
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
bool P1conform<Quadrilateral>::position_of_dof(int nrShapeFct, MathVector< RefDim >& value) const
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

namespace{

static const bool registered1 = ug::TrialSpaces::inst().register_trial_space(ug::TST_P1CONFORM, ug::P1conform<ug::Triangle>::inst());
static const bool registered2 = ug::TrialSpaces::inst().register_trial_space(ug::TST_P1CONFORM, ug::P1conform<ug::Quadrilateral>::inst());

}
