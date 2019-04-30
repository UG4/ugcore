/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_GRID__GRID_DIM_TRAITS__
#define __H__LIB_GRID__GRID_DIM_TRAITS__

#include <boost/mpl/list.hpp>
#include "lib_grid/grid_objects/grid_objects.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	Element types for different dimensions
////////////////////////////////////////////////////////////////////////////////

/**
 * This class template provides the basic element types and their topological
 * properties depending on the topological dimensionality. Besides that, it
 * provides the lists of the element types for every topological dimension.
 * While DimElemList returns the Element Types in the dimenion of the domain,
 * the list AllElemList returns all elements contained in the Domain-dimension
 * and the dimensions below.
 */
template <int dim> struct grid_dim_traits;

// 0d
template <> struct grid_dim_traits<0>
{
	typedef boost::mpl::list<RegularVertex> DimElemList;
	typedef boost::mpl::list<RegularVertex> AllElemList;

	typedef geometry_traits<Vertex>::const_iterator const_iterator;
	typedef geometry_traits<Vertex>::iterator iterator;

	typedef geometry_traits<Vertex>::grid_base_object grid_base_object;
	typedef geometry_traits<Vertex>::grid_base_object element_type;

	const static size_t MaxNumVerticesOfElem = 1;
};

// 1d
template <> struct grid_dim_traits<1>
{
	typedef boost::mpl::list<RegularEdge> DimElemList;
	typedef boost::mpl::list<RegularEdge> AllElemList;
	typedef boost::mpl::list<> ManifoldElemList;

	typedef geometry_traits<Edge>::const_iterator const_iterator;
	typedef geometry_traits<Edge>::iterator iterator;

	typedef geometry_traits<Edge>::grid_base_object grid_base_object;
	typedef geometry_traits<Edge>::grid_base_object element_type;
	typedef geometry_traits<Vertex>::grid_base_object side_type;

	const static size_t MaxNumVerticesOfElem = 2;
};

// 2d
template <> struct grid_dim_traits<2>
{
	typedef boost::mpl::list<Triangle, Quadrilateral> DimElemList;
	typedef boost::mpl::list<RegularEdge, Triangle, Quadrilateral> AllElemList;
	typedef boost::mpl::list<RegularEdge> ManifoldElemList;

	typedef geometry_traits<Face>::const_iterator const_iterator;
	typedef geometry_traits<Face>::iterator iterator;

	typedef geometry_traits<Face>::grid_base_object grid_base_object;
	typedef geometry_traits<Face>::grid_base_object element_type;
	typedef geometry_traits<Edge>::grid_base_object side_type;

	const static size_t MaxNumVerticesOfElem = 4;
};

// 3d
template <> struct grid_dim_traits<3>
{
	typedef boost::mpl::list<Tetrahedron, Prism, Pyramid, Hexahedron, Octahedron> DimElemList;
	typedef boost::mpl::list<RegularEdge, Triangle, Quadrilateral,
								Tetrahedron, Prism, Pyramid, Hexahedron, Octahedron> AllElemList;
	typedef boost::mpl::list<Triangle, Quadrilateral> ManifoldElemList;

	typedef geometry_traits<Volume>::const_iterator const_iterator;
	typedef geometry_traits<Volume>::iterator iterator;

	typedef geometry_traits<Volume>::grid_base_object grid_base_object;
	typedef geometry_traits<Volume>::grid_base_object element_type;
	typedef geometry_traits<Face>::grid_base_object side_type;

	const static size_t MaxNumVerticesOfElem = 8;
};

} // end namespace ug

#endif // __H__LIB_GRID__GRID_DIM_TRAITS__

/* End of File */
