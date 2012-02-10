// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 23.12.2011 (m,d,y)

#ifndef __H__UG__geometric_objects_0d__
#define __H__UG__geometric_objects_0d__

#include "../grid/grid.h"
#include "common/math/ugmath.h"
#include "common/assert.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	shared pipe sections vertex
///	These numbers define where in the vertex-section-container a vertex will be stored.
enum SharedPipeSectionVertex
{
	SPSVRT_NONE = -1,
	SPSVRT_VERTEX = 0,
	SPSVRT_HANGING_VERTEX = 1
};


////////////////////////////////////////////////////////////////////////
//	Vertex
///	A basic vertex-type
/**
 * This type represents the most commonly used vertex.
 *
 * \ingroup lib_grid_geometric_objects
 */
class UG_API Vertex : public VertexBase
{
	friend class Grid;
	public:
		inline static bool type_match(GeometricObject* pObj)	{return dynamic_cast<Vertex*>(pObj) != NULL;}

		virtual ~Vertex()	{}

		virtual GeometricObject* create_empty_instance() const	{return new Vertex;}

		virtual int shared_pipe_section() const	{return SPSVRT_VERTEX;}
		virtual int base_object_type_id() const	{return VERTEX;}
		virtual ReferenceObjectID reference_object_id() const {return ROID_VERTEX;}
};

template <>
class geometry_traits<Vertex>
{
	public:
		typedef GenericGeometricObjectIterator<Vertex*, VertexBaseIterator>			iterator;
		typedef ConstGenericGeometricObjectIterator<Vertex*, VertexBaseIterator,
														ConstVertexBaseIterator>	const_iterator;

		typedef VertexBase	geometric_base_object;

		enum
		{
			SHARED_PIPE_SECTION = SPSVRT_VERTEX,
			BASE_OBJECT_TYPE_ID = VERTEX
		};

		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_VERTEX;
};

typedef geometry_traits<Vertex>::iterator 		VertexIterator;
typedef geometry_traits<Vertex>::const_iterator	ConstVertexIterator;



////////////////////////////////////////////////////////////////////////
//	HangingVertex
///	A vertex appearing on edges or faces.
/**
 * HangingVertices appear during hangin-vertex-refinement.
 * They should be treated with care, since they are often referenced
 * by other elements, as i.e. \sa ConstrainingEdge.
 *
 * \ingroup lib_grid_geometric_objects
 */
class UG_API HangingVertex : public VertexBase
{
	friend class Grid;
	public:
		inline static bool type_match(GeometricObject* pObj)	{return dynamic_cast<HangingVertex*>(pObj) != NULL;}

		HangingVertex()	: m_pParent(NULL)	{}
		virtual ~HangingVertex()	{}

		virtual GeometricObject* create_empty_instance() const	{return new HangingVertex;}

		virtual int shared_pipe_section() const	{return SPSVRT_HANGING_VERTEX;}
		virtual int base_object_type_id() const	{return VERTEX;}
		virtual ReferenceObjectID reference_object_id() const {return ROID_VERTEX;}

		virtual bool is_constrained() const			{return true;}

		inline void set_parent(GeometricObject* pParent)	{m_pParent = pParent;}
		inline GeometricObject* get_parent()	{return m_pParent;}
		inline int get_parent_base_object_type_id()	{return m_pParent->base_object_type_id();}

		inline const vector2& get_local_coordinates() const	{return m_localCoord;}
		inline number get_local_coordinate_1() const		{return m_localCoord.x;}
		inline number get_local_coordinate_2() const		{return m_localCoord.y;}

		inline void set_local_coordinates(number x, number y)		{m_localCoord.x = x; m_localCoord.y = y;}
		inline void set_local_coordinates(const vector2& coords)	{m_localCoord = coords;}
		inline void set_local_coordinate_1(number coord) 	{m_localCoord.x = coord;}
		inline void set_local_coordinate_2(number coord) 	{m_localCoord.y = coord;}


	protected:
		GeometricObject*	m_pParent;
		vector2				m_localCoord;


};

template <>
class geometry_traits<HangingVertex>
{
	public:
		typedef GenericGeometricObjectIterator<HangingVertex*, VertexBaseIterator>			iterator;
		typedef ConstGenericGeometricObjectIterator<HangingVertex*, VertexBaseIterator,
																ConstVertexBaseIterator>	const_iterator;

		typedef VertexBase	geometric_base_object;

		enum
		{
			SHARED_PIPE_SECTION = SPSVRT_HANGING_VERTEX,
			BASE_OBJECT_TYPE_ID = VERTEX
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_VERTEX;
};

typedef geometry_traits<HangingVertex>::iterator 		HangingVertexIterator;
typedef geometry_traits<HangingVertex>::const_iterator	ConstHangingVertexIterator;


}//	end of namespace

#endif
