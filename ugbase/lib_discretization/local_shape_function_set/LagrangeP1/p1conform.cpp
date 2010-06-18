/*
 * p1conform.cpp
 *
 *  Created on: 16.02.2010
 *      Author: andreasvogel
 */


#include "p1conform.h"

namespace ug{

template<>
bool
P1conform<ReferenceTriangle>::
evaluate(int nrShapeFct, const position_type& locPos, shape_value_type& value) const
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
bool P1conform<ReferenceTriangle>::
evaluate_grad(int nrShapeFct, const position_type& locPos, grad_value_type& value) const
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
bool P1conform<ReferenceTriangle>::
position_of_dof(int nrShapeFct, position_type& value) const
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
bool
P1conform<ReferenceQuadrilateral>::
evaluate(int nrShapeFct, const position_type& locPos, shape_value_type& value) const
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
bool
P1conform<ReferenceQuadrilateral>::
evaluate_grad(int nrShapeFct, const position_type& locPos, grad_value_type& value) const
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
bool
P1conform<ReferenceQuadrilateral>::
position_of_dof(int nrShapeFct, position_type& value) const
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



template<>
bool
P1conform<ReferenceTetrahedron>::
evaluate(int nrShapeFct, const position_type& locPos, shape_value_type& value) const
{
	switch(nrShapeFct)
	{
		case 0: value = (1.-locPos.x-locPos.y-locPos.z); return true;
		case 1: value = locPos.x; return true;
		case 2: value = locPos.y; return true;
		case 3: value = locPos.z; return true;
		default: return false;
	}
	return true;
};

template<>
bool P1conform<ReferenceTetrahedron>::
evaluate_grad(int nrShapeFct, const position_type& locPos, grad_value_type& value) const
{
	switch(nrShapeFct)
	{
	case 0: value[0] = -1.0;
			value[1] = -1.0;
			value[2] = -1.0; return true;
	case 1: value[0] = 	1.0;
			value[1] =  0.0;
			value[2] =  0.0; return true;
	case 2: value[0] =  0.0;
			value[1] =  1.0;
			value[2] =  0.0; return true;
	case 3: value[0] =  0.0;
			value[1] =  0.0;
			value[2] =  1.0; return true;
	default: return false;
	}
	return true;
}

template<>
bool P1conform<ReferenceTetrahedron>::
position_of_dof(int nrShapeFct, position_type& value) const
{
	switch(nrShapeFct)
	{
	case 0: value[0] =  0.0;
			value[1] =  0.0;
			value[2] =  0.0; return true;
	case 1: value[0] = 	1.0;
			value[1] =  0.0;
			value[2] =  0.0; return true;
	case 2: value[0] =  0.0;
			value[1] =  1.0;
			value[2] =  0.0; return true;
	case 3: value[0] =  0.0;
			value[1] =  0.0;
			value[2] =  1.0; return true;
	default: return false;
	}
	return true;
}

}

