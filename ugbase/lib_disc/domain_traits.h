/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__DOMAIN_TRAITS__
#define __H__UG__LIB_DISC__DOMAIN_TRAITS__

#include <boost/mpl/for_each.hpp>
#include "lib_grid/grid_objects/grid_dim_traits.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	Element types for different dimensions
////////////////////////////////////////////////////////////////////////////////

/**
 * This Class provides boost::mpl::lists storing the type of elements used in
 * for the Domain. It can be used to control dimension dependent builds, where
 * not all template instantiations are available (e.g. a Hexahedron in 1d,2d, etc.)
 * While DimElemList returns the Element Types in the dimension of the domain,
 * the list AllElemList returns all elements contained in the Domain-dimension
 * and the dimensions below.
 */
template <int dim> struct domain_traits;

// 0d
template <> struct domain_traits<0> : grid_dim_traits<0> {};

// 1d
template <> struct domain_traits<1> : grid_dim_traits<1> {
    using position_type = MathVector<1>;
    using position_attachment_type = Attachment<position_type>;
    using position_accessor_type = Grid::VertexAttachmentAccessor<position_attachment_type>;
};

// 2d
template <> struct domain_traits<2> : grid_dim_traits<2> {
    using position_type = MathVector<2>;
    using position_attachment_type = Attachment<position_type>;
    using position_accessor_type = Grid::VertexAttachmentAccessor<position_attachment_type>;
};

// 3d
template <> struct domain_traits<3> : grid_dim_traits<3> {
    using position_type = MathVector<3>;
    using position_attachment_type = Attachment<position_type>;
    using position_accessor_type = Grid::VertexAttachmentAccessor<position_attachment_type>;
};


} // end namespace ug

#endif