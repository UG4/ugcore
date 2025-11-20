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

#ifndef __H__LIB_GRID__GEOMETRIC_OBJECT_COLLECTION__
#define __H__LIB_GRID__GEOMETRIC_OBJECT_COLLECTION__

#include <list>
#include "grid_base_objects.h"
#include "common/util/section_container.h"
#include "element_storage.h"
#include "grid_base_object_traits.h"

namespace ug
{

/// \addtogroup lib_grid_grid
/// @{

////////////////////////////////////////////////////////////////////////
//	GridObjectCollection
///	a helper class that holds a collection of possibly unconnected geometric-objects.
/**
 * This class is a simple helper class..
 * Its purpose is to make it easy to pass a collection of geometric-objects
 * to a function while maintaining the possibility to iterate over different
 * sub-types of geometric-objects seperatly.
 *
 * In contrary to \sa GridObjectCollection, the
 * GridObjectCollection allows access to the elements through
 * different levels.
 *
 * Please note that a GridObjectCollection is only valid as long as
 * the object from which you received the collection still exists.
 *
 * A GridObjectCollection can only be queried for iterators and 
 * element-counts. You may not insert new elements or remove
 * old ones (at least not directly).
 *
 * Classes that you can query for their GridObjectCollection
 * are for example ug::Grid, ug::MultiGrid, ug::SubsetHandler,
 * ug::MGSubsetHandler, ug::Selector, ug::MGSelector.
 *
 * As long as the object that provides the GridObjectCollection
 * is still valid, the GridObjectCollection will always hold the current
 * geometric objects of the source-object (i.e. a grid, a selector or a subset-handler),
 * as long as new objects are inserted into existing subsets (SubsetHandler) or
 * existing levels (MultiGrid). Insertion or removal of subsets or levels is
 * not reflected by the goc and can lead to severe errors.
 * Make sure to retrieve a new goc if such changes happened.
 *
 * Please note that a GridObjectCollection does not necessarily represent
 * a topological closed part of a grid.
 * A Collection can for example hold faces without their
 * associated vertices.
 *
 * How to use GridObjectCollection:
 * Once you retrieved an instance (let's call it goc) you can query it for
 * iterators like this:
 * VertexIterator iter = goc.vertices_begin(0);
 * or if you want to iterate over triangles of level 1 type the following:
 * TriangleIterator iter = goc.begin<Triangle>(1);
 *
 * if you want to get the number of hexahedrons in level 0 you would go like this:
 * uint numHexas = goc.num<Hexahedron>(0);
 */
class UG_API GridObjectCollection
{
	public:
	///	The traits class holds some important types for each element-type
		template <typename TElem>
		struct traits{
			using iterator = typename geometry_traits<TElem>::iterator;
			using const_iterator = typename geometry_traits<TElem>::const_iterator;
		};

	///	initializes the instance with an estimate of the number of levels.
	/**	The estimate does not have to match exactly. However, if it does
	 *  it makes things faster.*/
		GridObjectCollection(size_t levelEstimate = 1);
		
	///	initializes level 0 with the given sections.
		GridObjectCollection(ElementStorage<Vertex>::SectionContainer* vrtCon,
								ElementStorage<Edge>::SectionContainer* edgeCon,
								ElementStorage<Face>::SectionContainer* faceCon,
								ElementStorage<Volume>::SectionContainer* volCon);

	//	copy constructor.
		GridObjectCollection(const GridObjectCollection& mgoc);
		
		GridObjectCollection& operator = (const GridObjectCollection& mgoc);
		
	///	only used during creation by the methods that create the collection
		void add_level(ElementStorage<Vertex>::SectionContainer* vrtCon,
						ElementStorage<Edge>::SectionContainer* edgeCon,
						ElementStorage<Face>::SectionContainer* faceCon,
						ElementStorage<Volume>::SectionContainer* volCon);

	///	returns the number of levels
		inline size_t num_levels() const {return m_levels.size();}
		
	//	Iterators
	//	begin
	/**	returns the begin iterator for the specified level.
	 *	If no level is given iterators for level 0 are returned.*/
		template <typename TGeomObj>
		inline
		typename geometry_traits<TGeomObj>::iterator
		begin(size_t level = 0);

	//	end
	/**	returns the end iterator for the specified level.
	 *	If no level is given iterators for level 0 are returned.*/
		template <typename TGeomObj>
		inline
		typename geometry_traits<TGeomObj>::iterator
		end(size_t level = 0);

		inline VertexIterator	vertices_begin(size_t level = 0)	{return begin<Vertex>(level);}
		inline VertexIterator	vertices_end(size_t level = 0)		{return end<Vertex>(level);}
		inline EdgeIterator		edges_begin(size_t level = 0)		{return begin<Edge>(level);}
		inline EdgeIterator		edges_end(size_t level = 0)			{return end<Edge>(level);}
		inline FaceIterator			faces_begin(size_t level = 0)		{return begin<Face>(level);}
		inline FaceIterator			faces_end(size_t level = 0)			{return end<Face>(level);}
		inline VolumeIterator		volumes_begin(size_t level = 0)		{return begin<Volume>(level);}
		inline VolumeIterator		volumes_end(size_t level = 0)		{return end<Volume>(level);}

	//	const iterators
	//	begin
		template <typename TGeomObj>
		inline
		typename geometry_traits<TGeomObj>::const_iterator
		begin(size_t level = 0) const;

	//	end
		template <typename TGeomObj>
		inline
		typename geometry_traits<TGeomObj>::const_iterator
		end(size_t level = 0) const;

		inline ConstVertexIterator	vertices_begin(size_t level = 0) const	{return begin<Vertex>(level);}
		inline ConstVertexIterator	vertices_end(size_t level = 0) const	{return end<Vertex>(level);}
		inline ConstEdgeIterator	edges_begin(size_t level = 0) const		{return begin<Edge>(level);}
		inline ConstEdgeIterator	edges_end(size_t level = 0) const		{return end<Edge>(level);}
		inline ConstFaceIterator		faces_begin(size_t level = 0) const		{return begin<Face>(level);}
		inline ConstFaceIterator		faces_end(size_t level = 0) const		{return end<Face>(level);}
		inline ConstVolumeIterator		volumes_begin(size_t level = 0) const	{return begin<Volume>(level);}
		inline ConstVolumeIterator		volumes_end(size_t level = 0) const		{return end<Volume>(level);}
		
	//	element numbers
		template <typename TGeomObj>
		size_t num() const;
		
		inline size_t num_vertices() const	{return num<Vertex>();}
		inline size_t num_edges() const		{return num<Edge>();}
		inline size_t num_faces() const		{return num<Face>();}
		inline size_t num_volumes() const	{return num<Volume>();}
		
		template <typename TGeomObj>
		inline
		size_t num(size_t level) const;
		
		inline size_t num_vertices(size_t level) const	{return num<Vertex>(level);}
		inline size_t num_edges(size_t level) const		{return num<Edge>(level);}
		inline size_t num_faces(size_t level) const		{return num<Face>(level);}
		inline size_t num_volumes(size_t level) const	{return num<Volume>(level);}
		
	protected:
		void assign(const GridObjectCollection& goc);

		template <typename TGeomObj> inline
		const typename ElementStorage<typename geometry_traits<TGeomObj>::grid_base_object>::
		SectionContainer* get_container(size_t level) const;
		
		template <typename TGeomObj> inline
		typename ElementStorage<typename geometry_traits<TGeomObj>::grid_base_object>::
		SectionContainer* get_container(size_t level);
				
	protected:
		struct ContainerCollection{
			ContainerCollection() = default;
			ContainerCollection(ElementStorage<Vertex>::SectionContainer* vrtCon,
								ElementStorage<Edge>::SectionContainer* edgeCon,
								ElementStorage<Face>::SectionContainer* faceCon,
								ElementStorage<Volume>::SectionContainer* volCon);

			ElementStorage<Vertex>::SectionContainer* vrtContainer;
			ElementStorage<Edge>::SectionContainer* edgeContainer;
			ElementStorage<Face>::SectionContainer* faceContainer;
			ElementStorage<Volume>::SectionContainer* volContainer;
		};

		using ContainerVec = std::vector<ContainerCollection>;
		//using GOCVec = std::vector<GridObjectCollection>;

	protected:
		ContainerVec m_levels;
};

/// @}
}//	end of namespace

////////////////////////////////////////////////
//	include implementation
#include "grid_object_collection_impl.hpp"

#endif
