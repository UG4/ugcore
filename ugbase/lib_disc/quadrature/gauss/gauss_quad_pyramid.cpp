/*
 * Copyright (c) 2011-2013:  G-CSC, Goethe University Frankfurt
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
//  It provides the Gauss Quadratures for a reference pyramid.


#include "../quadrature.h"
#include "gauss_quad_pyramid.h"

namespace ug{

template <>
number GaussQuadBase<GaussQuadrature<ReferencePyramid, 2>, 3, 2, 8>::m_vWeight[8] =
{
1./3. * 0.12500000000000000000000000000000,
1./3. * 0.12500000000000000000000000000000,
1./3. * 0.12500000000000000000000000000000,
1./3. * 0.12500000000000000000000000000000,
1./3. * 0.12500000000000000000000000000000,
1./3. * 0.12500000000000000000000000000000,
1./3. * 0.12500000000000000000000000000000,
1./3. * 0.12500000000000000000000000000000
};

template <>
MathVector<3> GaussQuadBase<GaussQuadrature<ReferencePyramid, 2>, 3, 2, 8>::m_vPoint[8] =
{
MathVector<3>(0.58541020000000000000000000000000, 0.72360680000000000000000000000000, 0.13819660000000000000000000000000),
MathVector<3>(0.13819660000000000000000000000000, 0.72360680000000000000000000000000, 0.13819660000000000000000000000000),
MathVector<3>(0.13819660000000000000000000000000, 0.27639320000000000000000000000000, 0.58541020000000000000000000000000),
MathVector<3>(0.13819660000000000000000000000000, 0.27639320000000000000000000000000, 0.13819660000000000000000000000000),
MathVector<3>(0.72360680000000000000000000000000, 0.13819660000000000000000000000000, 0.13819660000000000000000000000000),
MathVector<3>(0.72360680000000000000000000000000, 0.58541020000000000000000000000000, 0.13819660000000000000000000000000),
MathVector<3>(0.27639320000000000000000000000000, 0.13819660000000000000000000000000, 0.58541020000000000000000000000000),
MathVector<3>(0.27639320000000000000000000000000, 0.13819660000000000000000000000000, 0.13819660000000000000000000000000)
};




template <>
FlexGaussQuadrature<ReferencePyramid>::FlexGaussQuadrature(int order)
{
	switch(order)
	{
	case 0:
	case 1:
	case 2:
		m_order = GaussQuadrature<ReferencePyramid, 2>::order();
		m_numPoints = GaussQuadrature<ReferencePyramid, 2>::size();
		m_pvPoint = GaussQuadrature<ReferencePyramid, 2>::points();
		m_pvWeight = GaussQuadrature<ReferencePyramid, 2>::weights();
		break;

	default: UG_THROW("Order "<<order<<" not available for GaussQuadrature of pyramid.");
	}
}
}; // namespace ug

