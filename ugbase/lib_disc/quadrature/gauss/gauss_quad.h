/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__QUADRATURE__GAUSS_QUAD__GAUSS_QUAD__
#define __H__UG__LIB_DISC__QUADRATURE__GAUSS_QUAD__GAUSS_QUAD__

#include "common/common.h"
#include "../quadrature.h"
#include "lib_disc/reference_element/reference_element.h"

namespace ug{

/// fixed order gauss quadrature
template <typename TRefElem, int order>
class GaussQuadrature;

/// wrapper to ease implementation
template <typename TImpl, int TDim, int TOrder, int TNip>
class GaussQuadBase
{
	public:
	/// Dimension of integration domain
		static constexpr size_t dim = TDim;

	/// Position Type in Reference Element Space
	using position_type = MathVector<dim>;

	///	Type of weights
	using weight_type = number;

	/// Order of quadrature rule
		static constexpr size_t p = TOrder;

	/// Number of integration points
		static constexpr size_t nip = TNip;

	public:
	/// number of integration points
		static size_t size() {return nip;}

	/// returns i'th integration point
		static const MathVector<dim>& point(size_t i)
			{UG_ASSERT(i < size(), "Wrong index"); return m_vPoint[i];}

	/// returns all positions in an array of size()
		static const MathVector<dim>* points() {return m_vPoint;}

	/// return the i'th weight
		static number weight(size_t i)
			{UG_ASSERT(i < size(), "Wrong index"); return m_vWeight[i];}

	/// returns all weights in an array of size()
		static const number* weights() {return m_vWeight;}

	/// returns the order
		static size_t order() {return p;}

	protected:
	/// integration points
		static MathVector<dim> m_vPoint[nip];

	/// weights
		static number m_vWeight[nip];
};

/// flexible order gauss quadrature
/**
 * Providing gauss quadrature for an reference element. This class wrapps a
 * static GaussQuadrature into the Quadrature interface.
 *
 * \tparam 		TRefElem		Reference Element Type
 */
template <typename TRefElem>
class FlexGaussQuadrature
	: public QuadratureRule<TRefElem::dim>
{
	public:
	///	Constructor
		FlexGaussQuadrature(int order);

	///	Destructor
		~FlexGaussQuadrature() {}
};

} // namespace ug

// include implementation
#include "gauss_quad_vertex.h"
#include "gauss_quad_edge.h"
#include "gauss_quad_triangle.h"
#include "gauss_quad_quadrilateral.h"
#include "gauss_quad_tetrahedron.h"
#include "gauss_quad_pyramid.h"
#include "gauss_quad_prism.h"
#include "gauss_quad_hexahedron.h"
#include "gauss_quad_octahedron.h"


#endif /* __H__UG__LIB_DISC__QUADRATURE__GAUSS_QUAD__GAUSS_QUAD__ */
