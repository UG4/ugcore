/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Dmitry Logashenko
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

/*
 * Classes for computing integral properties of geometric objects types in lists.
 */
#ifndef __H__UG__LIB_DISC__REFERENCE_ELEMENT__ELEMENT_LIST_TRAITS__
#define __H__UG__LIB_DISC__REFERENCE_ELEMENT__ELEMENT_LIST_TRAITS__

// further lib_disc headers
#include "lib_disc/reference_element/reference_element.h"

// Boost-C++ headers
#include <boost/mpl/transform_view.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/min_max.hpp>

namespace ug{

/// Class for computing integral properties of geometric objects types in lists
/**
 * This class gets the maximum number of corners, edges etc over all the geometric
 * elements in a given list. The computation is performed at the compile time
 * so that the result (e.g. element_list_traits<ElemList>::max_edges) is a static
 * constant.
 *
 * \tparam ElemList		list of geometric elements
 */
template <typename ElemList>
class element_list_traits
{
///	Metafunction class for counting corners in an element type
	struct mfc_num_corners_of_elem
	{
		template <typename TElem> class apply
		{
			using ref_elem_type = typename reference_element_traits<TElem>::reference_element_type;
			public:
				using type = boost::mpl::int_<ref_elem_type::numCorners>; ///< returned type (i.e. result of the metafunction)
		};
	};
	
///	Metafunction class for counting edges in an element type
	struct mfc_num_edges_of_elem
	{
		template <typename TElem> class apply
		{
			using ref_elem_type = typename reference_element_traits<TElem>::reference_element_type;
			public:
				using type = boost::mpl::int_<ref_elem_type::numEdges>; ///< returned type (i.e. result of the metafunction)
		};
	};
	
///	Metafunction class for counting faces in an element type
	struct mfc_num_faces_of_elem
	{
		template <typename TElem> class apply
		{
			using ref_elem_type = typename reference_element_traits<TElem>::reference_element_type;
			public:
				using type = boost::mpl::int_<ref_elem_type::numFaces>; ///< returned type (i.e. result of the metafunction)
		};
	};
	
///	Metafunction class for counting volumes in an element type
	struct mfc_num_volumes_of_elem
	{
		template <typename TElem> class apply
		{
			using ref_elem_type = typename reference_element_traits<TElem>::reference_element_type;
			public:
				using type = boost::mpl::int_<ref_elem_type::numVolumes>; ///< returned type (i.e. result of the metafunction)
		};
	};
	
///	Metafunction class for counting sides (i.e. edges or faces) in an element type
	struct mfc_num_sides_of_elem
	{
		template <typename TElem> class apply
		{
			using ref_elem_type = typename reference_element_traits<TElem>::reference_element_type;
			public:
				using type = boost::mpl::int_<ref_elem_type::numSides>; ///< returned type (i.e. result of the metafunction)
		};
	};
	
public:
	
/// Max. number of corners of the elements in the element list (as a constant)
	static constexpr int maxCorners = boost::mpl::fold
			< boost::mpl::transform_view<ElemList, mfc_num_corners_of_elem>,
				boost::mpl::int_<0>,
				boost::mpl::max<boost::mpl::_1,boost::mpl::_2>
			>::type::value;
	
/// Max. number of edges of the elements in the element list (as a constant)
	static constexpr int maxEdges = boost::mpl::fold
			< boost::mpl::transform_view<ElemList, mfc_num_edges_of_elem>,
				boost::mpl::int_<0>,
				boost::mpl::max<boost::mpl::_1,boost::mpl::_2>
			>::type::value;
	
/// Max. number of faces of the elements in the element list (as a constant)
	static constexpr int maxFaces = boost::mpl::fold
			< boost::mpl::transform_view<ElemList, mfc_num_faces_of_elem>,
				boost::mpl::int_<0>,
				boost::mpl::max<boost::mpl::_1,boost::mpl::_2>
			>::type::value;
	
/// Max. number of volumes of the elements in the element list (as a constant)
	static constexpr int maxVolumes = boost::mpl::fold
			< boost::mpl::transform_view<ElemList, mfc_num_volumes_of_elem>,
				boost::mpl::int_<0>,
				boost::mpl::max<boost::mpl::_1,boost::mpl::_2>
			>::type::value;
	
/// Max. number of sides (edges or faces) of the elements in the element list (as a constant)
	static constexpr int maxSides = boost::mpl::fold
			< boost::mpl::transform_view<ElemList, mfc_num_sides_of_elem>,
				boost::mpl::int_<0>,
				boost::mpl::max<boost::mpl::_1,boost::mpl::_2>
			>::type::value;
};

} // end namespace ug

#endif // __H__UG__LIB_DISC__REFERENCE_ELEMENT__ELEMENT_LIST_TRAITS__

/* End of File */
