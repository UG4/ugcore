/*
 * Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
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

//  This file is parsed from UG 3.9.
//  It provides the Gauss Quadratures for a reference vertex.

#ifndef __H__UG__LIB_DISC__QUADRATURE__GAUSS_QUAD__GAUSS_QUAD_VERTEX__
#define __H__UG__LIB_DISC__QUADRATURE__GAUSS_QUAD__GAUSS_QUAD_VERTEX__

#include "gauss_quad.h"

namespace ug{


// TODO: This class might be superfluous now that the specializations below exist.
class GaussQuadratureVertex : public QuadratureRule<0>
{
	public:
	/// Constructor
		GaussQuadratureVertex(){
			m_order = 20; // something large
			m_numPoints = 1;
			m_vPoint[0][0] = 0.000000000000000000;
			m_vWeight[0] =  1.000000000000000000;
			m_pvPoint = &m_vPoint[0];
			m_pvWeight = &m_vWeight[0];
		}

	protected:
	/// integration points
		MathVector<dim> m_vPoint[1];

	/// weights
		number m_vWeight[1];
};


template <>
class GaussQuadrature<ReferenceVertex, 0>
: public GaussQuadBase<GaussQuadrature<ReferenceVertex, 0>, 0, 0, 1> {};


template <>
class GaussQuadrature<ReferenceVertex, 1>
: public GaussQuadrature<ReferenceVertex, 0> {};

template <>
class GaussQuadrature<ReferenceVertex, 2>
: public GaussQuadrature<ReferenceVertex, 0> {};

template <>
class GaussQuadrature<ReferenceVertex, 3>
: public GaussQuadrature<ReferenceVertex, 0> {};

template <>
class GaussQuadrature<ReferenceVertex, 4>
: public GaussQuadrature<ReferenceVertex, 0> {};

template <>
class GaussQuadrature<ReferenceVertex, 5>
: public GaussQuadrature<ReferenceVertex, 0> {};

template <>
class GaussQuadrature<ReferenceVertex, 6>
: public GaussQuadrature<ReferenceVertex, 0> {};

template <>
class GaussQuadrature<ReferenceVertex, 7>
: public GaussQuadrature<ReferenceVertex, 0> {};

template <>
class GaussQuadrature<ReferenceVertex, 8>
: public GaussQuadrature<ReferenceVertex, 0> {};

template <>
class GaussQuadrature<ReferenceVertex, 9>
: public GaussQuadrature<ReferenceVertex, 0> {};

template <>
class GaussQuadrature<ReferenceVertex, 10>
: public GaussQuadrature<ReferenceVertex, 0> {};

template <>
class GaussQuadrature<ReferenceVertex, 11>
: public GaussQuadrature<ReferenceVertex, 0> {};

template <>
class GaussQuadrature<ReferenceVertex, 12>
: public GaussQuadrature<ReferenceVertex, 0> {};

template <>
class GaussQuadrature<ReferenceVertex, 13>
: public GaussQuadrature<ReferenceVertex, 0> {};

template <>
class GaussQuadrature<ReferenceVertex, 14>
: public GaussQuadrature<ReferenceVertex, 0> {};

template <>
class GaussQuadrature<ReferenceVertex, 15>
: public GaussQuadrature<ReferenceVertex, 0> {};

template <>
class GaussQuadrature<ReferenceVertex, 16>
: public GaussQuadrature<ReferenceVertex, 0> {};

template <>
class GaussQuadrature<ReferenceVertex, 17>
: public GaussQuadrature<ReferenceVertex, 0> {};

template <>
class GaussQuadrature<ReferenceVertex, 18>
: public GaussQuadrature<ReferenceVertex, 0> {};

template <>
class GaussQuadrature<ReferenceVertex, 19>
: public GaussQuadrature<ReferenceVertex, 0> {};

}; // namespace ug


#endif