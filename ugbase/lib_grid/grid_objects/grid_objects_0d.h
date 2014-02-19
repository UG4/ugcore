// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 23.12.2011 (m,d,y)

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
class UG_API RegularVertex : public VertexBase
{
	friend class Grid;
	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<RegularVertex*>(pObj) != NULL;}

		virtual ~RegularVertex()	{}

		virtual GridObject* create_empty_instance() const	{return new RegularVertex;}

		virtual int container_section() const	{return CSVRT_REGULAR_VERTEX;}
		virtual ReferenceObjectID reference_object_id() const {return ROID_VERTEX;}
};

template <>
class geometry_traits<RegularVertex>
{
	public:
		typedef GenericGridObjectIterator<RegularVertex*, VertexBaseIterator>			iterator;
		typedef ConstGenericGridObjectIterator<RegularVertex*, VertexBaseIterator,
														ConstVertexBaseIterator>	const_iterator;

		typedef VertexBase	grid_base_object;

		enum
		{
			CONTAINER_SECTION = CSVRT_REGULAR_VERTEX,
			BASE_OBJECT_ID = VERTEX
		};

		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_VERTEX;
};

typedef geometry_traits<RegularVertex>::iterator 		RegularVertexIterator;
typedef geometry_traits<RegularVertex>::const_iterator	ConstRegularVertexIterator;



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
class UG_API ConstrainedVertex : public VertexBase
{
	friend class Grid;
	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<ConstrainedVertex*>(pObj) != NULL;}

		ConstrainedVertex()	: m_constrainingObj(NULL), m_parentBaseObjectId(-1)	{}
		virtual ~ConstrainedVertex()	{}

		virtual GridObject* create_empty_instance() const	{return new ConstrainedVertex;}

		virtual int container_section() const	{return CSVRT_CONSTRAINED_VERTEX;}
		virtual ReferenceObjectID reference_object_id() const {return ROID_VERTEX;}

		virtual bool is_constrained() const			{return true;}

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
		vector2				m_localCoord;
		int					m_parentBaseObjectId;


};

template <>
class geometry_traits<ConstrainedVertex>
{
	public:
		typedef GenericGridObjectIterator<ConstrainedVertex*, VertexBaseIterator>			iterator;
		typedef ConstGenericGridObjectIterator<ConstrainedVertex*, VertexBaseIterator,
																ConstVertexBaseIterator>	const_iterator;

		typedef VertexBase	grid_base_object;

		enum
		{
			CONTAINER_SECTION = CSVRT_CONSTRAINED_VERTEX,
			BASE_OBJECT_ID = VERTEX
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_VERTEX;
};

typedef geometry_traits<ConstrainedVertex>::iterator 		ConstrainedVertexIterator;
typedef geometry_traits<ConstrainedVertex>::const_iterator	ConstConstrainedVertexIterator;


}//	end of namespace

#endif
