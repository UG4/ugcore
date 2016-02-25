/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */


#include "./lagrangep1.h"

namespace ug{

/// \cond HIDDEN_SYMBOLS

///////////////////////////////////////
// ReferenceVertex
///////////////////////////////////////

template<>
LagrangeP1<ReferenceVertex>::shape_type
LagrangeP1<ReferenceVertex>::
shape(size_t i, const MathVector<dim>& x) const
{
	return 1.0;
};

template<>
void
LagrangeP1<ReferenceVertex>::
grad(grad_type& value, size_t i, const MathVector<dim>& x) const
{
}

template<>
bool LagrangeP1<ReferenceVertex>::
position(size_t i, MathVector<dim>& value) const
{
	return true;
}

///////////////////////////////////////
// ReferenceEdge
///////////////////////////////////////

template<>
LagrangeP1<ReferenceEdge>::shape_type
LagrangeP1<ReferenceEdge>::
shape(size_t i, const MathVector<dim>& x) const
{
	switch(i)
	{
	case 0: return (1.-x[0]);
	case 1: return x[0];
	default: UG_THROW("LagrangeP1: Invalid shape fct index: "<<i);
	}
};

template<>
void
LagrangeP1<ReferenceEdge>::
grad(grad_type& value, size_t i, const MathVector<dim>& x) const
{
	switch(i)
	{
	case 0: value[0] = -1.0; break;
	case 1: value[0] = 	1.0; break;
	default: UG_THROW("LagrangeP1: Invalid shape fct index: "<<i);
	}
}

template<>
bool LagrangeP1<ReferenceEdge>::
position(size_t i, MathVector<dim>& value) const
{
	switch(i)
	{
	case 0: value[0] =  0.0; return true;
	case 1: value[0] = 	1.0; return true;
	default: UG_THROW("LagrangeP1: Invalid shape fct index: "<<i);
	}
	return true;
}

///////////////////////////////////////
// ReferenceTriangle
///////////////////////////////////////

template<>
LagrangeP1<ReferenceTriangle>::shape_type
LagrangeP1<ReferenceTriangle>::
shape(size_t i, const MathVector<dim>& x) const
{
	switch (i)
	{
	case 0: return(1.0-x[0]-x[1]);
	case 1: return(x[0]);
	case 2: return(x[1]);
	default: UG_THROW("LagrangeP1: Invalid shape fct index: "<<i);
	}
};

template<>
void
LagrangeP1<ReferenceTriangle>::
grad(grad_type& value, size_t i, const MathVector<dim>& x) const
{
	switch(i)
	{
	case 0: value[0] = -1.0;
			value[1] = -1.0; break;
	case 1: value[0] = 	1.0;
			value[1] =  0.0; break;
	case 2: value[0] =  0.0;
			value[1] =  1.0; break;
	default: UG_THROW("LagrangeP1: Invalid shape fct index: "<<i);
	}
}

template<>
bool LagrangeP1<ReferenceTriangle>::
position(size_t i, MathVector<dim>& value) const
{
	switch(i)
	{
	case 0: value[0] =  0.0;
			value[1] =  0.0; return true;
	case 1: value[0] = 	1.0;
			value[1] =  0.0; return true;
	case 2: value[0] =  0.0;
			value[1] =  1.0; return true;
	default: UG_THROW("LagrangeP1: Invalid shape fct index: "<<i);
	}
	return true;
}

///////////////////////////////////////
// ReferenceQuadrilateral
///////////////////////////////////////

template<>
LagrangeP1<ReferenceQuadrilateral>::shape_type
LagrangeP1<ReferenceQuadrilateral>::
shape(size_t i, const MathVector<dim>& x) const
{
	switch(i)
	{
	case 0: return((1.0-x[0])*(1.0-x[1]));
	case 1: return(x[0]*(1.0-x[1]));
	case 2: return(x[0]*x[1]);
	case 3: return((1.0-x[0])*x[1]);
	default: UG_THROW("LagrangeP1: Invalid shape fct index: "<<i);
	}
};

template<>
void
LagrangeP1<ReferenceQuadrilateral>::
grad(grad_type& value, size_t i, const MathVector<dim>& x) const
{
	switch(i)
	{
	case 0: value[0] = -1.0 + x[1];
			value[1] = -1.0 + x[0]; break;
	case 1: value[0] = 	1.0 - x[1];
			value[1] =      - x[0]; break;
	case 2: value[0] =        x[1];
			value[1] =        x[0]; break;
	case 3: value[0] =      - x[1];
			value[1] =  1.0 - x[0]; break;
	default: UG_THROW("LagrangeP1: Invalid shape fct index: "<<i);
	}
}

template<>
bool
LagrangeP1<ReferenceQuadrilateral>::
position(size_t i, MathVector<dim>& value) const
{
	switch(i)
	{
	case 0: value[0] =  0.0;
			value[1] =  0.0; return true;
	case 1: value[0] = 	1.0;
			value[1] =  0.0; return true;
	case 2: value[0] =  1.0;
			value[1] =  1.0; return true;
	case 3: value[0] =  0.0;
			value[1] =  1.0; return true;
	default: UG_THROW("LagrangeP1: Invalid shape fct index: "<<i);
	}
	return true;
}

///////////////////////////////////////
// ReferenceTetrahedron
///////////////////////////////////////

template<>
LagrangeP1<ReferenceTetrahedron>::shape_type
LagrangeP1<ReferenceTetrahedron>::
shape(size_t i, const MathVector<dim>& x) const
{
	switch(i)
	{
	case 0: return(1.0-x[0]-x[1]-x[2]);
	case 1: return(x[0]);
	case 2: return(x[1]);
	case 3: return(x[2]);
	default: UG_THROW("LagrangeP1: Invalid shape fct index: "<<i);
	}
};

template<>
void
LagrangeP1<ReferenceTetrahedron>::
grad(grad_type& value, size_t i, const MathVector<dim>& x) const
{
	switch(i)
	{
	case 0: value[0] = -1.0;
			value[1] = -1.0;
			value[2] = -1.0; break;
	case 1: value[0] = 	1.0;
			value[1] =  0.0;
			value[2] =  0.0; break;
	case 2: value[0] =  0.0;
			value[1] =  1.0;
			value[2] =  0.0; break;
	case 3: value[0] =  0.0;
			value[1] =  0.0;
			value[2] =  1.0; break;
	default: UG_THROW("LagrangeP1: Invalid shape fct index: "<<i);
	}
}

template<>
bool LagrangeP1<ReferenceTetrahedron>::
position(size_t i, MathVector<dim>& value) const
{
	switch(i)
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
	default: UG_THROW("LagrangeP1: Invalid shape fct index: "<<i);
	}
	return true;
}

///////////////////////////////////////
// ReferencePyramid
///////////////////////////////////////

template<>
LagrangeP1<ReferencePyramid>::shape_type
LagrangeP1<ReferencePyramid>::
shape(size_t i, const MathVector<dim>& x) const
{
	switch(i)
	{
	  case 0 :
		if (x[0] > x[1])
		  return((1.0-x[0])*(1.0-x[1]) + x[2]*(x[1]-1.0));
		else
		  return((1.0-x[0])*(1.0-x[1]) + x[2]*(x[0]-1.0));
	  case 1 :
		if (x[0] > x[1])
		  return(x[0]*(1.0-x[1])       - x[2]*x[1]);
		else
		  return(x[0]*(1.0-x[1])       - x[2]*x[0]);
	  case 2 :
		if (x[0] > x[1])
		  return(x[0]*x[1]             + x[2]*x[1]);
		else
		  return(x[0]*x[1]             + x[2]*x[0]);
	  case 3 :
		if (x[0] > x[1])
		  return((1.0-x[0])*x[1]       - x[2]*x[1]);
		else
		  return((1.0-x[0])*x[1]       - x[2]*x[0]);
	  case 4 : return(x[2]);
	default: UG_THROW("LagrangeP1: Invalid shape fct index: "<<i);
	}
};

template<>
void
LagrangeP1<ReferencePyramid>::
grad(grad_type& value, size_t i, const MathVector<dim>& x) const
{
	switch(i)
	{
	  case 0:
		if (x[0] > x[1])
		  {
			value[0] = -(1.0-x[1]);
			value[1] = -(1.0-x[0]) + x[2];
			value[2] = -(1.0-x[1]);
			break;
		  }
		else
		  {
			value[0] = -(1.0-x[1]) + x[2];
			value[1] = -(1.0-x[0]);
			value[2] = -(1.0-x[0]);
			break;
		  }
	  case 1:
		if (x[0] > x[1])
		  {
			value[0] = (1.0-x[1]);
			value[1] = -x[0] - x[2];
			value[2] = -x[1];
			break;
		  }
		else
		  {
			value[0] = (1.0-x[1]) - x[2];
			value[1] = -x[0];
			value[2] = -x[0];
			break;
		  }
	  case 2:
		if (x[0] > x[1])
		  {
			value[0] = x[1];
			value[1] = x[0] + x[2];
			value[2] = x[1];
			break;
		  }
		else
		  {
			value[0] = x[1] + x[2];
			value[1] = x[0];
			value[2] = x[0];
			break;
		  }
	  case 3:
		if (x[0] > x[1])
		  {
			value[0] = -x[1];
			value[1] = 1.0-x[0] - x[2];
			value[2] = -x[1];
			break;
		  }
		else
		  {
			value[0] = -x[1] - x[2];
			value[1] = 1.0-x[0];
			value[2] = -x[0];
			break;
		  }
      case 4:
		value[0] = 0.0;
		value[1] = 0.0;
		value[2] = 1.0;
		break;
	default: UG_THROW("LagrangeP1: Invalid shape fct index: "<<i);
	}
}

template<>
bool LagrangeP1<ReferencePyramid>::
position(size_t i, MathVector<dim>& value) const
{
	static const DimReferenceElement<3>& refElem
		= ReferenceElementProvider::get<3>(ROID_PYRAMID);

	value = refElem.corner(i);
	return true;
}

///////////////////////////////////////
// ReferencePrism
///////////////////////////////////////

template<>
LagrangeP1<ReferencePrism>::shape_type
LagrangeP1<ReferencePrism>::
shape(size_t i, const MathVector<dim>& x) const
{
	switch(i)
	{
	case 0 : return((1.0-x[0]-x[1])*(1.0-x[2]));
	case 1 : return(x[0]*(1.0-x[2]));
	case 2 : return(x[1]*(1.0-x[2]));
	case 3 : return((1.0-x[0]-x[1])*x[2]);
	case 4 : return(x[0]*x[2]);
	case 5 : return(x[1]*x[2]);
	default: UG_THROW("LagrangeP1: Invalid shape fct index: "<<i);
	}
};

template<>
void
LagrangeP1<ReferencePrism>::
grad(grad_type& value, size_t i, const MathVector<dim>& x) const
{
	switch(i)
	{
	  case 0:
		value[0] = -(1.0-x[2]);
		value[1] = -(1.0-x[2]);
		value[2] = -(1.0-x[0]-x[1]);
		break;
	  case 1:
		value[0] = (1.0-x[2]);
		value[1] = 0.0;
		value[2] = -x[0];
		break;
	  case 2:
		value[0] = 0.0;
		value[1] = (1.0-x[2]);
		value[2] = -x[1];
		break;
	  case 3:
		value[0] = -x[2];
		value[1] = -x[2];
		value[2] = 1.0-x[0]-x[1];
      break;
    case 4:
		value[0] = x[2];
		value[1] = 0.0;
		value[2] = x[0];
		break;
    case 5:
		value[0] = 0.0;
		value[1] = x[2];
		value[2] = x[1];
		break;
	default: UG_THROW("LagrangeP1: Invalid shape fct index: "<<i);
	}
}

template<>
bool LagrangeP1<ReferencePrism>::
position(size_t i, MathVector<dim>& value) const
{
	static const DimReferenceElement<3>& refElem
		= ReferenceElementProvider::get<3>(ROID_PRISM);

	value = refElem.corner(i);
	return true;
}

///////////////////////////////////////
// ReferenceHexahedron
///////////////////////////////////////

template<>
LagrangeP1<ReferenceHexahedron>::shape_type
LagrangeP1<ReferenceHexahedron>::
shape(size_t i, const MathVector<dim>& x) const
{
	switch(i)
	{
	case 0: return((1.0-x[0])*(1.0-x[1])*(1.0-x[2]));
	case 1: return((x[0])*(1.0-x[1])*(1.0-x[2]));
	case 2: return((x[0])*(x[1])*(1.0-x[2]));
	case 3: return((1.0-x[0])*(x[1])*(1.0-x[2]));
	case 4: return((1.0-x[0])*(1.0-x[1])*(x[2]));
	case 5: return((x[0])*(1.0-x[1])*(x[2]));
	case 6: return((x[0])*(x[1])*(x[2]));
	case 7: return((1.0-x[0])*(x[1])*(x[2]));
	default: UG_THROW("LagrangeP1: Invalid shape fct index: "<<i);
	}
};

template<>
void
LagrangeP1<ReferenceHexahedron>::
grad(grad_type& value, size_t i, const MathVector<dim>& x) const
{
	switch(i)
	{
	  case 0:
		value[0] = -(1.0-x[1])*(1.0-x[2]);
		value[1] = -(1.0-x[0])*(1.0-x[2]);
		value[2] = -(1.0-x[0])*(1.0-x[1]);
		break;
	  case 1:
		value[0] = (1.0-x[1])*(1.0-x[2]);
		value[1] = -(x[0])*(1.0-x[2]);
		value[2] = -(x[0])*(1.0-x[1]);
		break;
	  case 2:
		value[0] = (x[1])*(1.0-x[2]);
		value[1] = (x[0])*(1.0-x[2]);
		value[2] = -x[0]*x[1];
		break;
	  case 3:
		value[0] = -(x[1])*(1.0-x[2]);
		value[1] = (1.0-x[0])*(1.0-x[2]);
		value[2] = -(1.0-x[0])*(x[1]);
      break;
    case 4:
		value[0] = -(1.0-x[1])*(x[2]);
		value[1] = -(1.0-x[0])*(x[2]);
		value[2] = (1.0-x[0])*(1.0-x[1]);
		break;
	  case 5:
		value[0] = (1.0-x[1])*x[2];
		value[1] = -(x[0])*x[2];
		value[2] = (x[0])*(1.0-x[1]);
		break;
	  case 6:
		value[0] = (x[1])*x[2];
		value[1] = (x[0])*x[2];
		value[2] = x[0]*x[1];
		break;
	  case 7:
		value[0] = -(x[1])*x[2];
		value[1] = (1.0-x[0])*x[2];
		value[2] = (1.0-x[0])*x[1];
      break;
	default: UG_THROW("LagrangeP1: Invalid shape fct index: "<<i);
	}
}

template<>
bool LagrangeP1<ReferenceHexahedron>::
position(size_t i, MathVector<dim>& value) const
{
	static const DimReferenceElement<3>& refElem
		= ReferenceElementProvider::get<3>(ROID_HEXAHEDRON);

	value = refElem.corner(i);
	return true;
}

///////////////////////////////////////
// ReferenceOctahedron
///////////////////////////////////////

template<>
LagrangeP1<ReferenceOctahedron>::shape_type
LagrangeP1<ReferenceOctahedron>::
shape(size_t i, const MathVector<dim>& x) const
{
//	the octahedral shape functions correspond to the tetrahedral ones,
//	but locally piecewise linear on a subdivision of the octahedron
//	into 4 sub-tetrahedrons.
	const number z_sgn 	= (x[2] < 0) ? -1.0 : 1.0;

	switch(i)
	{
	  case 0 :
		if (x[2] < 0)
		  return(-1.0*x[2]);
		else
		  return(0.0);
	  case 1 :
		if (x[0] > x[1])
		  return(1.0-x[0]-z_sgn*x[2]);
		else
		  return(1.0-x[1]-z_sgn*x[2]);
	  case 2 :
		if (x[0] > x[1])
		  return(x[0]-x[1]);
		else
	      return(0.0);
	  case 3 :
		if (x[0] > x[1])
		  return(x[1]);
		else
		  return(x[0]);
	  case 4 :
		if (x[0] > x[1])
		  return(0.0);
		else
		  return(x[1]-x[0]);
	  case 5 :
		if (x[2] < 0)
		  return(0.0);
		else
		  return(x[2]);
	default: UG_THROW("LagrangeP1: Invalid shape fct index: "<<i);
	}
};

template<>
void
LagrangeP1<ReferenceOctahedron>::
grad(grad_type& value, size_t i, const MathVector<dim>& x) const
{
//	the octahedral shape functions correspond to the tetrahedral ones,
//	but locally piecewise linear on a subdivision of the octahedron
//	into 4 sub-tetrahedrons.
	const number z_sgn 	= (x[2] < 0) ? -1.0 : 1.0;

	switch(i)
	{
	  case 0:
		if (x[2] < 0.0)
		{
			value[0] = 0.0;
			value[1] = 0.0;
			value[2] = -1.0;
			break;
		}
		else
		{
			value[0] = 0.0;
			value[1] = 0.0;
			value[2] = 0.0;
			break;
		}
	  case 1:
		if (x[0] > x[1])
		  {
			value[0] = -1.0;
			value[1] =  0.0;
			value[2] = -z_sgn;
			break;
		  }
		else
		  {
			value[0] =  0.0;
			value[1] = -1.0;
			value[2] = -z_sgn;
			break;
		  }
	  case 2:
		if (x[0] > x[1])
		  {
			value[0] =  1.0;
			value[1] = -1.0;
			value[2] =  0.0;
			break;
		  }
		else
		  {
			value[0] = 0.0;
			value[1] = 0.0;
			value[2] = 0.0;
			break;
		  }
	  case 3:
		if (x[0] > x[1])
		  {
			value[0] =  0.0;
			value[1] =  1.0;
			value[2] =  0.0;
			break;
		  }
		else
		  {
			value[0] =  1.0;
			value[1] =  0.0;
			value[2] =  0.0;
			break;
		  }
	  case 4:
		if (x[0] > x[1])
		  {
			value[0] =  0.0;
			value[1] =  0.0;
			value[2] =  0.0;
			break;
		  }
		else
		  {
			value[0] = -1.0;
			value[1] =  1.0;
			value[2] =  0.0;
			break;
		  }
	  case 5:
		if (x[2] < 0.0)
		{
			value[0] = 0.0;
			value[1] = 0.0;
			value[2] = 0.0;
			break;
		}
		else
		{
			value[0] = 0.0;
			value[1] = 0.0;
			value[2] = 1.0;
			break;
		}
	default: UG_THROW("LagrangeP1: Invalid shape fct index: "<<i);
	}
}

template<>
bool LagrangeP1<ReferenceOctahedron>::
position(size_t i, MathVector<dim>& value) const
{
	static const DimReferenceElement<3>& refElem
		= ReferenceElementProvider::get<3>(ROID_OCTAHEDRON);

	value = refElem.corner(i);
	return true;
}

/// \endcond

}

