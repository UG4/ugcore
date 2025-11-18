/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__grid_objects_1d__
#define __H__UG__grid_objects_1d__

#include "../grid/grid.h"
#include "common/math/ugmath.h"
#include "common/assert.h"
#include "grid_objects_0d.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
///	These numbers define where in the edge-section-container an edge will be stored.
/**	The order of the constants must not be changed! Algorithms may exist that rely on it.*/
enum EdgeContainerSections
{
	CSEDGE_NONE = -1,
	CSEDGE_REGULAR_EDGE = 0,
	CSEDGE_CONSTRAINED_EDGE = 1,
	CSEDGE_CONSTRAINING_EDGE = 2
};



////////////////////////////////////////////////////////////////////////
//	RegularEdge
///	Edges connect two vertices.
/**
 * The most commonly used edge-type.
 *
 * \ingroup lib_grid_grid_objects
 */
class UG_API RegularEdge : public Edge
{
	friend class Grid;
	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<RegularEdge*>(pObj) != nullptr;}

		RegularEdge() = default;
		RegularEdge(Vertex* v1, Vertex* v2)
			{
				m_vertices[0] = v1;
				m_vertices[1] = v2;
			}

		RegularEdge(const EdgeDescriptor& descriptor)
			{
				m_vertices[0] = descriptor.vertex(0);
				m_vertices[1] = descriptor.vertex(1);
			}

		virtual ~RegularEdge() = default;

		virtual GridObject* create_empty_instance() const	{return new RegularEdge;}

		virtual int container_section() const	{return CSEDGE_REGULAR_EDGE;}
		virtual ReferenceObjectID reference_object_id() const {return ROID_EDGE;}

	///	virtual refine. Returns pointers to Edge.
	/**
	 * create 2 new edges, connecting the original edges end-points with vrtNew.
	 * \sa Edge::refine.
	 */
		virtual bool refine(std::vector<Edge*>& vNewEdgesOut,
							Vertex* newVertex,
							Vertex** pSubstituteVrts = nullptr);

//TODO:	Think about this method. It is not safe!
	///	non virtual refine. Returns pointers to RegularEdge.
	/**
	 * create 2 new edges, connecting the original edges end-points with vrtNew.
	 * \sa Edge::refine.
	 */
		bool refine(std::vector<RegularEdge*>& vNewEdgesOut,
					Vertex* newVertex,
					Vertex** pSubstituteVrts = nullptr);
};

template <>
class geometry_traits<RegularEdge>
{
	public:
		using iterator = GenericGridObjectIterator<RegularEdge*, EdgeIterator>;
		using const_iterator = ConstGenericGridObjectIterator<RegularEdge*, EdgeIterator,
			ConstEdgeIterator>;

		using Descriptor = EdgeDescriptor;
		using grid_base_object = Edge;

		enum
		{
			CONTAINER_SECTION = CSEDGE_REGULAR_EDGE,
			BASE_OBJECT_ID = EDGE
		};
		static constexpr ReferenceObjectID REFERENCE_OBJECT_ID = ROID_EDGE;
};

using RegularEdgeIterator = geometry_traits<RegularEdge>::iterator;
using ConstRegularEdgeIterator = geometry_traits<RegularEdge>::const_iterator;



////////////////////////////////////////////////////////////////////////
//	ConstrainedEdge
///	This edge is a sub-edge of a \sa ConstrainingEdge.
/**
 * Edges of this type appear during hanging-vertex-refinement.
 * Treat them with care, since they are referenced by other objects,
 * i.e. \sa ConstrainingEdge
 *
 * \ingroup lib_grid_grid_objects
 */
class UG_API ConstrainedEdge : public Edge
{
	friend class Grid;
	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<ConstrainedEdge*>(pObj) != nullptr;}

		ConstrainedEdge() : m_pConstrainingObject(nullptr), m_parentBaseObjectId(-1)	{}
		ConstrainedEdge(Vertex* v1, Vertex* v2) :
			m_pConstrainingObject(nullptr),
			m_parentBaseObjectId(-1)
			{
				m_vertices[0] = v1;
				m_vertices[1] = v2;
			}

		ConstrainedEdge(const EdgeDescriptor& descriptor) :
			m_pConstrainingObject(nullptr),
			m_parentBaseObjectId(-1)
			{
				m_vertices[0] = descriptor.vertex(0);
				m_vertices[1] = descriptor.vertex(1);
			}

		virtual ~ConstrainedEdge()
		{
			if(m_pConstrainingObject)
				m_pConstrainingObject->remove_constraint_link(this);
		}

		virtual GridObject* create_empty_instance() const	{return new ConstrainedEdge;}

		virtual int container_section() const	{return CSEDGE_CONSTRAINED_EDGE;}
		virtual ReferenceObjectID reference_object_id() const {return ROID_EDGE;}

		virtual bool is_constrained() const			{return true;}

		virtual void remove_constraint_link(const Edge* e)
		{
			if(m_pConstrainingObject == static_cast<const GridObject*>(e)){
				m_pConstrainingObject = nullptr;
			}
		}

		virtual void remove_constraint_link(const Face* f)
		{
			if(m_pConstrainingObject == static_cast<const GridObject*>(f)){
				m_pConstrainingObject = nullptr;
			}
		}

	///	virtual refine. Returns pointers to Edge.
	/**
	 * create 2 new constrained edges, connecting the original edges end-points with vrtNew.
	 * please note that the caller is responsible to resolve any conflicts arising with
	 * existing links of constrained/constraining objects.
	 * \sa Edge::refine.
	 */
		virtual bool refine(std::vector<Edge*>& vNewEdgesOut,
							Vertex* newVertex,
							Vertex** pSubstituteVrts = nullptr);

//TODO:	Think about this method. It is not safe!
	///	non virtual refine. Returns pointers to ConstrainedEdge.
	/**
	 * create 2 new constrained edges, connecting the original edges end-points with vrtNew.
	 * please note that the caller is responsible to resolve any conflicts arising with
	 * existing links of constrained/constraining objects.
	 * \sa Edge::refine
	 */
		bool refine(std::vector<ConstrainedEdge*>& vNewEdgesOut,
					Vertex* newVertex,
					Vertex** pSubstituteVrts = nullptr);

		inline void set_constraining_object(GridObject* pObj)
		{
			m_pConstrainingObject = pObj;
			if(pObj)
				m_parentBaseObjectId = pObj->base_object_id();
		}

		inline GridObject* get_constraining_object()	{return m_pConstrainingObject;}

		inline int get_parent_base_object_id()				{return m_parentBaseObjectId;}
		inline void set_parent_base_object_id(int id)
		{
			if((m_parentBaseObjectId != -1) && (m_parentBaseObjectId != id)){
				UG_THROW("Bad parent base object id specified! The given id"
						" has to match the id of the constraining object if that"
						" is present. Call this method only, if no constraining"
						" object has been set!");
			}
			m_parentBaseObjectId = id;
		}


	protected:
		GridObject*	m_pConstrainingObject;
		int			m_parentBaseObjectId;
};

template <>
class geometry_traits<ConstrainedEdge>
{
	public:
	using iterator = GenericGridObjectIterator<ConstrainedEdge*, EdgeIterator>;
		using const_iterator = ConstGenericGridObjectIterator<ConstrainedEdge*, EdgeIterator,
			ConstEdgeIterator>;

		using Descriptor = EdgeDescriptor;
		using grid_base_object = Edge;

		enum
		{
			CONTAINER_SECTION = CSEDGE_CONSTRAINED_EDGE,
			BASE_OBJECT_ID = EDGE
		};
		static constexpr ReferenceObjectID REFERENCE_OBJECT_ID = ROID_EDGE;
};

using ConstrainedEdgeIterator = geometry_traits<ConstrainedEdge>::iterator;
using ConstConstrainedEdgeIterator = geometry_traits<ConstrainedEdge>::const_iterator;



////////////////////////////////////////////////////////////////////////
//	ConstrainingEdge
///	contains elements of type \sa HangingVertex and \sa ConstrainedEdge
/**
 * Edges of this type appear during hanging-vertex-refinement.
 * Treat with care.
 *
 * \ingroup lib_grid_grid_objects
 */
class UG_API ConstrainingEdge : public Edge
{
	friend class Grid;
	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<ConstrainingEdge*>(pObj) != nullptr;}

		ConstrainingEdge()	{}
		ConstrainingEdge(Vertex* v1, Vertex* v2)
			{
				m_vertices[0] = v1;
				m_vertices[1] = v2;
				m_constrainedVertices.reserve(1);
				m_constrainedEdges.reserve(2);
			}

		ConstrainingEdge(const EdgeDescriptor& descriptor)
			{
				m_vertices[0] = descriptor.vertex(0);
				m_vertices[1] = descriptor.vertex(1);
				m_constrainedVertices.reserve(1);
				m_constrainedEdges.reserve(2);
			}

		virtual ~ConstrainingEdge()
		{
			for(size_t i = 0; i < m_constrainedVertices.size(); ++i){
				m_constrainedVertices[i]->remove_constraint_link(this);
			}

			for(size_t i = 0; i < m_constrainedEdges.size(); ++i){
				m_constrainedEdges[i]->remove_constraint_link(this);
			}
		}

		virtual GridObject* create_empty_instance() const	{return new ConstrainingEdge;}

		virtual int container_section() const	{return CSEDGE_CONSTRAINING_EDGE;}
		virtual ReferenceObjectID reference_object_id() const {return ROID_EDGE;}

		virtual bool is_constraining() const	{return true;}

		
		virtual void remove_constraint_link(const Vertex* vrt)
		{
			unconstrain_object(vrt);
		}

		virtual void remove_constraint_link(const Edge* e)
		{
			unconstrain_object(e);
		}


	///	virtual refine. Returns pointers to Edge.
	/**
	 * create 2 new constraining edges, connecting the original edges end-points with vrtNew.
	 * The refine has no effect on constrained objects. The user has to manually copy them.
	 * \sa Edge::refine.
	 */
		virtual bool refine(std::vector<Edge*>& vNewEdgesOut,
							Vertex* newVertex,
							Vertex** pSubstituteVrts = nullptr);

//TODO:	Think about this method. It is not safe!
	///	non virtual refine. Returns pointers to ConstrainingEdge.
	/**
	 * create 2 new constraining edges, connecting the original edges end-points with vrtNew.
	 * The refine has no effect on constrained objects. The user has to manually copy them.
	 * \sa Edge::refine.
	 */
		bool refine(std::vector<ConstrainingEdge*>& vNewEdgesOut,
						Vertex* newVertex,
						Vertex** pSubstituteVrts = nullptr);


		inline void add_constrained_object(Vertex* pObj)
			{
				UG_ASSERT(!is_constrained_object(pObj), "vertex is already constrained by this edge.");
					m_constrainedVertices.push_back(pObj);
			}

		inline void add_constrained_object(Edge* pObj)
			{
				UG_ASSERT(!is_constrained_object(pObj), "edge is already constrained by this edge.");
					m_constrainedEdges.push_back(pObj);
			}

		inline bool is_constrained_object(Vertex* vrt)
			{
				std::vector<Vertex*>::iterator iter = find(m_constrainedVertices.begin(),
																m_constrainedVertices.end(), vrt);
				return iter != m_constrainedVertices.end();
			}

		inline bool is_constrained_object(Edge* edge)
			{
				std::vector<Edge*>::iterator iter = find(m_constrainedEdges.begin(),
																m_constrainedEdges.end(), edge);
				return iter != m_constrainedEdges.end();
			}

		inline void unconstrain_object(const Vertex* vrt)
			{
				std::vector<Vertex*>::iterator iter = find(m_constrainedVertices.begin(),
																m_constrainedVertices.end(), vrt);
				if(iter != m_constrainedVertices.end())
					m_constrainedVertices.erase(iter);
			}

		inline void unconstrain_object(const Edge* edge)
			{
				std::vector<Edge*>::iterator iter = find(m_constrainedEdges.begin(),
																m_constrainedEdges.end(), edge);
				if(iter != m_constrainedEdges.end())
					m_constrainedEdges.erase(iter);
			}

		inline void clear_constrained_vertices()	{m_constrainedVertices.clear();}
		inline void clear_constrained_edges()		{m_constrainedEdges.clear();}
		inline void clear_constrained_objects()
			{
				clear_constrained_vertices();
				clear_constrained_edges();
			}

	//	NUMBER OF CONSTRAINED ELEMENTS
		inline size_t num_constrained_vertices() const	{return m_constrainedVertices.size();}
		inline size_t num_constrained_edges() const		{return m_constrainedEdges.size();}

		template <class TElem> size_t num_constrained() const;


	//	ACCESS TO CONSTRAINED ELEMENTS
		Vertex* constrained_vertex(size_t ind) const
		{
			UG_ASSERT(ind < m_constrainedVertices.size(), "bad index");
			return m_constrainedVertices[ind];
		}

		Edge* constrained_edge(size_t ind) const
		{
			UG_ASSERT(ind < m_constrainedEdges.size(), "bad index");
			return m_constrainedEdges[ind];
		}

		template <class TElem> TElem* constrained(size_t ind) const;

	protected:
		std::vector<Vertex*>	m_constrainedVertices;
		std::vector<Edge*>		m_constrainedEdges;
};

template <>
class geometry_traits<ConstrainingEdge>
{
	public:
	using iterator = GenericGridObjectIterator<ConstrainingEdge*, EdgeIterator>;
		using const_iterator = ConstGenericGridObjectIterator<ConstrainingEdge*, EdgeIterator,
			ConstEdgeIterator>;

		using Descriptor = EdgeDescriptor;
		using grid_base_object = Edge;

		enum
		{
			CONTAINER_SECTION = CSEDGE_CONSTRAINING_EDGE,
			BASE_OBJECT_ID = EDGE
		};
		static constexpr ReferenceObjectID REFERENCE_OBJECT_ID = ROID_EDGE;
};

using ConstrainingEdgeIterator = geometry_traits<ConstrainingEdge>::iterator;
using ConstConstrainingEdgeIterator = geometry_traits<ConstrainingEdge>::const_iterator;

}//	end of namespace

#endif
