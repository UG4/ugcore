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

#ifndef __H__UG__grid_objects_0d__
#define __H__UG__grid_objects_0d__

#include "../grid/grid.h"
#include "common/math/ugmath.h"
#include "common/assert.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
///	These numbers define where in the vertex-section-container a vertex will be stored.
/**	The order of the constants must not be changed! Algorithms may exist that rely on it.*/
enum VertexContainerSections
{
	CSVRT_NONE = -1,
	CSVRT_REGULAR_VERTEX = 0,
	CSVRT_CONSTRAINED_VERTEX = 1
};


////////////////////////////////////////////////////////////////////////
//	RegularVertex
///	A basic vertex-type
/**
 * This type represents the most commonly used vertex.
 *
 * \ingroup lib_grid_grid_objects
 */
class UG_API RegularVertex : public Vertex
{
	friend class Grid;
	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<RegularVertex*>(pObj) != nullptr;}

		virtual ~RegularVertex()	{}

		virtual GridObject* create_empty_instance() const	{return new RegularVertex;}

		virtual int container_section() const	{return CSVRT_REGULAR_VERTEX;}
		virtual ReferenceObjectID reference_object_id() const {return ROID_VERTEX;}
};

template <>
class geometry_traits<RegularVertex>
{
	public:
		using iterator = GenericGridObjectIterator<RegularVertex*, VertexIterator>;
		using const_iterator = ConstGenericGridObjectIterator<RegularVertex*, VertexIterator,
			ConstVertexIterator>;

		using grid_base_object = Vertex;

		enum
		{
			CONTAINER_SECTION = CSVRT_REGULAR_VERTEX,
			BASE_OBJECT_ID = VERTEX
		};

		static constexpr ReferenceObjectID REFERENCE_OBJECT_ID = ROID_VERTEX;
};

using RegularVertexIterator = geometry_traits<RegularVertex>::iterator;
using ConstRegularVertexIterator = geometry_traits<RegularVertex>::const_iterator;



////////////////////////////////////////////////////////////////////////
//	ConstrainedVertex
///	A vertex appearing on edges or faces.
/**
 * ConstrainedVertices appear during hanging-vertex-refinement.
 * They should be treated with care, since they are often referenced
 * by other elements, as e.g. by a ConstrainingEdge.
 *
 * \ingroup lib_grid_grid_objects
 */
class UG_API ConstrainedVertex : public Vertex
{
	friend class Grid;
	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<ConstrainedVertex*>(pObj) != nullptr;}

		ConstrainedVertex()	: m_constrainingObj(nullptr), m_parentBaseObjectId(-1)	{}
		virtual ~ConstrainedVertex()
		{
			if(m_constrainingObj)
				m_constrainingObj->remove_constraint_link(this);
		}

		virtual GridObject* create_empty_instance() const	{return new ConstrainedVertex;}

		virtual int container_section() const	{return CSVRT_CONSTRAINED_VERTEX;}
		virtual ReferenceObjectID reference_object_id() const {return ROID_VERTEX;}

		virtual bool is_constrained() const			{return true;}

		virtual void remove_constraint_link(const Edge* e)
		{
			if(m_constrainingObj == static_cast<const GridObject*>(e)){
				m_constrainingObj = nullptr;
			}
		}

		virtual void remove_constraint_link(const Face* f)
		{
			if(m_constrainingObj == static_cast<const GridObject*>(f)){
				m_constrainingObj = nullptr;
			}
		}

		inline void set_constraining_object(GridObject* constrObj)
		{
			m_constrainingObj = constrObj;
			if(constrObj)
				m_parentBaseObjectId = constrObj->base_object_id();
		}

		inline GridObject* get_constraining_object()	{return m_constrainingObj;}
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

		inline const vector2& get_local_coordinates() const	{return m_localCoord;}
		inline number get_local_coordinate_1() const		{return m_localCoord.x();}
		inline number get_local_coordinate_2() const		{return m_localCoord.y();}

		inline void set_local_coordinates(number x, number y)		{m_localCoord.x() = x; m_localCoord.y() = y;}
		inline void set_local_coordinates(const vector2& coords)	{m_localCoord = coords;}
		inline void set_local_coordinate_1(number coord) 			{m_localCoord.x() = coord;}
		inline void set_local_coordinate_2(number coord) 			{m_localCoord.y() = coord;}

	protected:
		GridObject*	m_constrainingObj;
		vector2		m_localCoord;
		int			m_parentBaseObjectId;
};

template <>
class geometry_traits<ConstrainedVertex>
{
	public:
		using iterator = GenericGridObjectIterator<ConstrainedVertex*, VertexIterator>;
		using const_iterator = ConstGenericGridObjectIterator<ConstrainedVertex*, VertexIterator,
			ConstVertexIterator>;

		using grid_base_object = Vertex;

		enum
		{
			CONTAINER_SECTION = CSVRT_CONSTRAINED_VERTEX,
			BASE_OBJECT_ID = VERTEX
		};
		static constexpr ReferenceObjectID REFERENCE_OBJECT_ID = ROID_VERTEX;
};

using ConstrainedVertexIterator = geometry_traits<ConstrainedVertex>::iterator;
using ConstConstrainedVertexIterator = geometry_traits<ConstrainedVertex>::const_iterator;


}//	end of namespace

#endif
