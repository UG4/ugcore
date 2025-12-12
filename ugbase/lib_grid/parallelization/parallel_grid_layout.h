/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__LIB_GRID__PARALLEL_GRID_LAYOUT__
#define __H__LIB_GRID__PARALLEL_GRID_LAYOUT__

//#include <vector>
#include <list>
#include <map>
//#include <algorithm>
#include "pcl/pcl_communication_structs.h"
#include "lib_grid/grid/grid_base_objects.h"


//	specialize pcl::type_traits for Vertex, Edge, Face and Volume
namespace pcl {
///	Vertex interfaces and layouts store elements of type Vertex*
template <>
struct type_traits<ug::Vertex>
{
	using Elem = ug::Vertex*;
};

///	Edge interfaces and layouts store elements of type Vertex*
template <>
struct type_traits<ug::Edge>
{
	using Elem = ug::Edge*;
};

///	Face interfaces and layouts store elements of type Vertex*
template <>
struct type_traits<ug::Face>
{
	using Elem = ug::Face*;
};

///	Volume interfaces and layouts store elements of type Vertex*
template <>
struct type_traits<ug::Volume>
{
	using Elem = ug::Volume*;
};

}//	end of namespace pcl

namespace ug {

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///	The types of interface-entries.
/**	INT_H_MASTER and INT_H_SLAVE describe (horizontal) connections between
 *	nodes on one level in a grid-hierarchy. They are used to communicate
 *	data between neighbours.
 *
 *	INT_V_MASTER and INT_V_SLAVE describe (vertical) connections
 *	between nodes on different levels of a grid. They are used to
 *	communicate data between parents and children. They are only used
 *	for multigrids.
 *
 *	Note that the type parameter in DistributionInterfaceEntry is currently
 *	restricted to 4 bytes only!
 *
 *	developer note: Introducing INT_HORIZONTAL, INT_VERTICAL, INT_MASTER,
 *					INT_SLAVE and building or combinations of those
 *					would not make sense! Think about INT_H_MASTER | INT_V_SLAVE
 *					in this case (it would not be clear whether the master
 *					was in the H or in the V interface).
 */
enum InterfaceNodeTypes : byte_t
{
	INT_NONE =	0,
	INT_H_MASTER = 1,		///< horizontal master node
	INT_H_SLAVE = 1<<1,		///< horizontal slave node
	INT_V_MASTER = 1<<2,	///< vertical master node
	INT_V_SLAVE = 1<<3,		///< vertical slave node
};

//	declare vertex-, edge-, face- and volume-layouts
//	we're using std::list as interface-element container, since we
//	require interface-element-iterators that stay valid even if the
//	interface is altered.
//	Make sure that those layouts match the ones in GridLayoutMap.
using VertexLayout = pcl::MultiLevelLayout< pcl::OrderedInterface<Vertex, std::list> >;
using EdgeLayout = pcl::MultiLevelLayout< pcl::OrderedInterface<Edge, std::list> >;
using FaceLayout = pcl::MultiLevelLayout< pcl::OrderedInterface<Face, std::list> >;
using VolumeLayout = pcl::MultiLevelLayout< pcl::OrderedInterface<Volume, std::list> >;


////////////////////////////////////////////////////////////////////////
//	GridLayoutMap
///	lets you access layouts by type and key
/**
 * The GridLayoutMap helps you to organize your layouts
 * (e.g. master- and slave-layouts).
 *
 * You may query layouts for Vertex, Edge, Face and Volume.
 *
 * You may use a LayoutMap as follows:
 *
 * \code
 * GridLayoutMap layoutMap;
 * assert(!layoutMap.has_layout<Vertex>(0));
 * VertexLayout& layout = layoutMap.get_layout<Vertex>(0);
 * assert(layoutMap.has_layout<Vertex>(0));
 * \endcode
 *
 * To get associated types you may use the GridLayoutMap::Types array:
 * \code
 * GridLayoutMap::Types<Vertex>::Layout l = layoutMap.get_layout<Vertex>(0);
 * \endcode
 *
 * The Types struct is very useful when it comes to using a LayoutMap in
 * template code, too.
 */
class GridLayoutMap
{
	public:
		using Key = int;

	///	defines the types that are used by a LayoutMap for a given TType.
		template <typename TType>
		struct Types
		{
			using Interface = pcl::OrderedInterface<TType, std::list>;
			using Layout = pcl::MultiLevelLayout<Interface>;
			using Element = typename Interface::Element;
			using Map = std::map<Key, Layout>;
		};

	public:
	///	checks whether the layout associated with the given key exists for the given type.
		template <typename TType>
		[[nodiscard]] bool
		has_layout(const Key& key) const;

	///	creates the required layout if it doesn't exist already.
		template <typename TType>
		typename Types<TType>::Layout&
		get_layout(const Key& key);

		template <typename TType>
		const typename Types<TType>::Layout&
		get_layout(const Key& key) const;

	///	begin-iterator to the layout-map for the given type.
	/**	iter.first will return the key, iter.second the layout
	 *	(of type LayoutMap::Types<TType>::Layout).
	 *	\{ */
		template <typename TType>
		typename Types<TType>::Map::iterator
		layouts_begin();

		template <typename TType>
		typename Types<TType>::Map::const_iterator
		layouts_begin() const;
	/** \} */

	///	end-iterator to the layout-map for the given type.
	/**	iter.first will return the key, iter.second the layout
	 *	(of type LayoutMap::Types<TType>::Layout).
	 *	\{ */
		template <typename TType>
		typename Types<TType>::Map::iterator
		layouts_end();

		template <typename TType>
		typename Types<TType>::Map::const_iterator
		layouts_end() const;
	/** \} */

	///	erases the specified layout
	/**	returns an iterator to the next layout.*/
		template <typename TType>
		typename Types<TType>::Map::iterator
		erase_layout(typename Types<TType>::Map::iterator iter);
								
	///	erases the specified layout if it exists
		template <typename TType>
		void erase_layout(const Key& key);

		void clear();

	///	removes empty interfaces.
		void remove_empty_interfaces();

	private:
		template <typename TType>
		inline typename Types<TType>::Map&
		get_layout_map();
		
		template <typename TType>
		inline const typename Types<TType>::Map&
		get_layout_map() const;

	///	the argument is only a dummy to allow to choose the right method at compile time
	// \{
		inline Types<Vertex>::Map&
		get_layout_map(Vertex*);

		inline const Types<Vertex>::Map&
		get_layout_map(Vertex*) const;

		inline Types<Edge>::Map&
		get_layout_map(Edge*);

		inline const Types<Edge>::Map&
		get_layout_map(Edge*) const;

		inline Types<Face>::Map&
		get_layout_map(Face*);

		inline const Types<Face>::Map&
		get_layout_map(Face*) const;

		inline Types<Volume>::Map&
		get_layout_map(Volume*);

		inline const Types<Volume>::Map&
		get_layout_map(Volume*) const;
	// \}
	
	private:
		Types<Vertex>::Map m_vertexLayoutMap;
		Types<Edge>::Map m_edgeLayoutMap;
		Types<Face>::Map m_faceLayoutMap;
		Types<Volume>::Map m_volumeLayoutMap;
};

/// @}
}//	end of namespace

////////////////////////////////
//	include implementation
#include "parallel_grid_layout_impl.hpp"

#endif
