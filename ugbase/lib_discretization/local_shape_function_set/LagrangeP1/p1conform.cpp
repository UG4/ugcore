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
P1conform<ReferenceEdge>::
evaluate(int nrShapeFct, const position_type& locPos, shape_value_type& value) const
{
	switch(nrShapeFct)
	{
		case 0: value = (1.-locPos[0]); return true;
		case 1: value = locPos[0]; return true;
		default: return false;
	}
	return true;
};

template<>
bool P1conform<ReferenceEdge>::
evaluate_grad(int nrShapeFct, const position_type& locPos, grad_value_type& value) const
{
	switch(nrShapeFct)
	{
	case 0: value[0] = -1.0; return true;
	case 1: value[0] = 	1.0; return true;
	default: return false;
	}
	return true;
}

template<>
bool P1conform<ReferenceEdge>::
position_of_dof(int nrShapeFct, position_type& value) const
{
	switch(nrShapeFct)
	{
	case 0: value[0] =  0.0; return true;
	case 1: value[0] = 	1.0; return true;
	default: return false;
	}
	return true;
}


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

template<>
bool
P1conform<ReferencePyramid>::
evaluate(int nrShapeFct, const position_type& locPos, shape_value_type& value) const
{
	switch(nrShapeFct)
	{
	  case 0 :
		if (locPos.x > locPos.y)
		{ value = ((1.0-locPos.x)*(1.0-locPos.y) - locPos.z*(1.0-locPos.y)); return true;}
		else
		{ value = ((1.0-locPos.x)*(1.0-locPos.y) - locPos.z*(1.0-locPos.x)); return true;}
	  case 1 :
		if (locPos.x > locPos.y)
		{ value = (locPos.x*(1.0-locPos.y)       - locPos.z*locPos.y); return true;}
		else
		{ value = (locPos.x*(1.0-locPos.y)       - locPos.z*locPos.x); return true;}
	  case 2 :
		if (locPos.x > locPos.y)
		{ value = (locPos.x*locPos.y             + locPos.z*locPos.y); return true;}
		else
		{ value = (locPos.x*locPos.y             + locPos.z*locPos.x); return true;}
	  case 3 :
		if (locPos.x > locPos.y)
		{ value = ((1.0-locPos.x)*locPos.y       - locPos.z*locPos.y); return true;}
		else
		{ value = ((1.0-locPos.x)*locPos.y       - locPos.z*locPos.x); return true;}
	  case 4 : { value = locPos.z; return true; }
	  default: return false;
	}
	return true;
};

template<>
bool P1conform<ReferencePyramid>::
evaluate_grad(int nrShapeFct, const position_type& locPos, grad_value_type& value) const
{
	switch(nrShapeFct)
	{
	  case 0:
		if (locPos.x > locPos.y)
		  {
			value[0] = -(1.0-locPos.y);
			value[1] = -(1.0-locPos.x) + locPos.z;
			value[2] = -(1.0-locPos.y);
			return true;
		  }
		else
		  {
			value[0] = -(1.0-locPos.y) + locPos.z;
			value[1] = -(1.0-locPos.x);
			value[2] = -(1.0-locPos.x);
			return true;
		  }
	  case 1:
		if (locPos.x > locPos.y)
		  {
			value[0] = (1.0-locPos.y);
			value[1] = -locPos.x - locPos.z;
			value[2] = -locPos.y;
			return true;
		  }
		else
		  {
			value[0] = (1.0-locPos.y) - locPos.z;
			value[1] = -locPos.x;
			value[2] = -locPos.x;
			return true;
		  }
	  case 2:
		if (locPos.x > locPos.y)
		  {
			value[0] = locPos.y;
			value[1] = locPos.x + locPos.z;
			value[2] = locPos.y;
			return true;
		  }
		else
		  {
			value[0] = locPos.y + locPos.z;
			value[1] = locPos.x;
			value[2] = locPos.x;
			return true;
		  }
	  case 3:
		if (locPos.x > locPos.y)
		  {
			value[0] = -locPos.y;
			value[1] = 1.0-locPos.x - locPos.z;
			value[2] = -locPos.y;
			return true;
		  }
		else
		  {
			value[0] = -locPos.y - locPos.z;
			value[1] = 1.0-locPos.x;
			value[2] = -locPos.x;
			return true;
		  }
      case 4:
		value[0] = 0.0;
		value[1] = 0.0;
		value[2] = 1.0;
		return true;
	default: return false;
	}
	return true;
}

template<>
bool P1conform<ReferencePyramid>::
position_of_dof(int nrShapeFct, position_type& value) const
{
	static const DimReferenceElement<3>& refElem
		= DimReferenceElementFactory<3>::get_reference_element(ROID_PYRAMID);

	value = refElem.corner(nrShapeFct);
	return true;
}


template<>
bool
P1conform<ReferencePrism>::
evaluate(int nrShapeFct, const position_type& locPos, shape_value_type& value) const
{
	switch(nrShapeFct)
	{
	case 0: value = (1.0-locPos.x-locPos.y)*(1.0-locPos.z); return true;
	case 1: value = locPos.x*(1.0-locPos.z); return true;
	case 2: value = locPos.y*(1.0-locPos.z); return true;
	case 3: value = (1.0-locPos.x-locPos.y)*locPos.z; return true;
	case 4: value = locPos.x*locPos.z; return true;
	case 5: value = locPos.y*locPos.z; return true;
	default: return false;
	}
	return true;
};

template<>
bool P1conform<ReferencePrism>::
evaluate_grad(int nrShapeFct, const position_type& locPos, grad_value_type& value) const
{
	switch(nrShapeFct)
	{
	  case 0:
		value[0] = -(1.0-locPos.z);
		value[1] = -(1.0-locPos.z);
		value[2] = -(1.0-locPos.x-locPos.y);
		return true;
	  case 1:
		value[0] = (1.0-locPos.z);
		value[1] = 0.0;
		value[2] = -locPos.x;
		return true;
	  case 2:
		value[0] = 0.0;
		value[1] = (1.0-locPos.z);
		value[2] = -locPos.y;
		return true;
	  case 3:
		value[0] = -locPos.z;
		value[1] = -locPos.z;
		value[2] = 1.0-locPos.x-locPos.y;
      return true;
    case 4:
		value[0] = locPos.z;
		value[1] = 0.0;
		value[2] = locPos.x;
		return true;
    case 5:
		value[0] = 0.0;
		value[1] = locPos.z;
		value[2] = locPos.y;
		return true;
	default: return false;
	}
	return true;
}

template<>
bool P1conform<ReferencePrism>::
position_of_dof(int nrShapeFct, position_type& value) const
{
	static const DimReferenceElement<3>& refElem
		= DimReferenceElementFactory<3>::get_reference_element(ROID_PRISM);

	value = refElem.corner(nrShapeFct);
	return true;
}



template<>
bool
P1conform<ReferenceHexahedron>::
evaluate(int nrShapeFct, const position_type& locPos, shape_value_type& value) const
{
	switch(nrShapeFct)
	{
	  case 0 : value = ((1.0-locPos.x)*(1.0-locPos.y)*(1.0-locPos.z)); return true;
	  case 1 : value = (locPos.x*(1.0-locPos.y)*(1.0-locPos.z)); return true;
	  case 2 : value = (locPos.x*locPos.y*(1.0-locPos.z)); return true;
	  case 3 : value = ((1.0-locPos.x)*locPos.y*(1.0-locPos.z)); return true;
	  case 4 : value = ((1.0-locPos.x)*(1.0-locPos.y)*locPos.z); return true;
	  case 5 : value = (locPos.x*(1.0-locPos.y)*locPos.z); return true;
	  case 6 : value = (locPos.x*locPos.y*locPos.z); return true;
	  case 7 : value = ((1.0-locPos.x)*locPos.y*locPos.z); return true;
	default: return false;
	}
	return true;
};

template<>
bool P1conform<ReferenceHexahedron>::
evaluate_grad(int nrShapeFct, const position_type& locPos, grad_value_type& value) const
{
	switch(nrShapeFct)
	{
	  case 0:
		value[0] = -(1.0-locPos.y)*(1.0-locPos.z);
		value[1] = -(1.0-locPos.x)*(1.0-locPos.z);
		value[2] = -(1.0-locPos.x)*(1.0-locPos.y);
		return true;
	  case 1:
		value[0] = (1.0-locPos.y)*(1.0-locPos.z);
		value[1] = -(locPos.x)*(1.0-locPos.z);
		value[2] = -(locPos.x)*(1.0-locPos.y);
		return true;
	  case 2:
		value[0] = (locPos.y)*(1.0-locPos.z);
		value[1] = (locPos.x)*(1.0-locPos.z);
		value[2] = -locPos.x*locPos.y;
		return true;
	  case 3:
		value[0] = -(locPos.y)*(1.0-locPos.z);
		value[1] = (1.0-locPos.x)*(1.0-locPos.z);
		value[2] = -(1.0-locPos.x)*(locPos.y);
      return true;
    case 4:
		value[0] = -(1.0-locPos.y)*(locPos.z);
		value[1] = -(1.0-locPos.x)*(locPos.z);
		value[2] = (1.0-locPos.x)*(1.0-locPos.y);
		return true;
	  case 5:
		value[0] = (1.0-locPos.y)*locPos.z;
		value[1] = -(locPos.x)*locPos.z;
		value[2] = (locPos.x)*(1.0-locPos.y);
		return true;
	  case 6:
		value[0] = (locPos.y)*locPos.z;
		value[1] = (locPos.x)*locPos.z;
		value[2] = locPos.x*locPos.y;
		return true;
	  case 7:
		value[0] = -(locPos.y)*locPos.z;
		value[1] = (1.0-locPos.x)*locPos.z;
		value[2] = (1.0-locPos.x)*locPos.y;
      return true;
	default: return false;
	}
	return true;
}

template<>
bool P1conform<ReferenceHexahedron>::
position_of_dof(int nrShapeFct, position_type& value) const
{
	static const DimReferenceElement<3>& refElem
		= DimReferenceElementFactory<3>::get_reference_element(ROID_HEXAHEDRON);

	value = refElem.corner(nrShapeFct);
	return true;
}


}

