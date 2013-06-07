// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 23.12.2011 (m,d,y)

#ifndef __H__UG__geometric_objects_1d__
#define __H__UG__geometric_objects_1d__

#include "../grid/grid.h"
#include "common/math/ugmath.h"
#include "common/assert.h"
#include "geometric_objects_0d.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
///	These numbers define where in the edge-section-container an edge will be stored.
/**	The order of the constants must not be changed! Algorithms may exist that rely on it.*/
enum EdgeContainerSections
{
	CSEDGE_NONE = -1,
	CSEDGE_EDGE = 0,
	CSEDGE_CONSTRAINED_EDGE = 1,
	CSEDGE_CONSTRAINING_EDGE = 2
};



////////////////////////////////////////////////////////////////////////
//	Edge
///	Edges connect two vertices.
/**
 * The most commonly used edge-type.
 *
 * \ingroup lib_grid_geometric_objects
 */
class UG_API Edge : public EdgeBase
{
	friend class Grid;
	public:
		inline static bool type_match(GeometricObject* pObj)	{return dynamic_cast<Edge*>(pObj) != NULL;}

		Edge()	{}
		Edge(VertexBase* v1, VertexBase* v2)
			{
				m_vertices[0] = v1;
				m_vertices[1] = v2;
			}

		Edge(const EdgeDescriptor& descriptor)
			{
				m_vertices[0] = descriptor.vertex(0);
				m_vertices[1] = descriptor.vertex(1);
			}

		virtual ~Edge()	{}

		virtual GeometricObject* create_empty_instance() const	{return new Edge;}

		virtual int container_section() const	{return CSEDGE_EDGE;}
		virtual ReferenceObjectID reference_object_id() const {return ROID_EDGE;}

	///	virtual refine. Returns pointers to EdgeBase.
	/**
	 * create 2 new edges, connecting the original edges end-points with vrtNew.
	 * \sa EdgeBase::refine.
	 */
		virtual bool refine(std::vector<EdgeBase*>& vNewEdgesOut,
							VertexBase* newVertex,
							VertexBase** pSubstituteVrts = NULL);

//TODO:	Think about this method. It is not safe!
	///	non virtual refine. Returns pointers to Edge.
	/**
	 * create 2 new edges, connecting the original edges end-points with vrtNew.
	 * \sa EdgeBase::refine.
	 */
		bool refine(std::vector<Edge*>& vNewEdgesOut,
					VertexBase* newVertex,
					VertexBase** pSubstituteVrts = NULL);
};

template <>
class geometry_traits<Edge>
{
	public:
		typedef GenericGeometricObjectIterator<Edge*, EdgeBaseIterator>			iterator;
		typedef ConstGenericGeometricObjectIterator<Edge*, EdgeBaseIterator,
														ConstEdgeBaseIterator>	const_iterator;

		typedef EdgeDescriptor	Descriptor;
		typedef EdgeBase		geometric_base_object;

		enum
		{
			CONTAINER_SECTION = CSEDGE_EDGE,
			BASE_OBJECT_ID = EDGE
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_EDGE;
};

typedef geometry_traits<Edge>::iterator 		EdgeIterator;
typedef geometry_traits<Edge>::const_iterator 	ConstEdgeIterator;



////////////////////////////////////////////////////////////////////////
//	ConstrainedEdge
///	This edge is a sub-edge of a \sa ConstrainingEdge.
/**
 * Edges of this type appear during hanging-vertex-refinement.
 * Treat them with care, since they are referenced by other objects,
 * i.e. \sa ConstrainingEdge
 *
 * \ingroup lib_grid_geometric_objects
 */
class UG_API ConstrainedEdge : public EdgeBase
{
	friend class Grid;
	public:
		inline static bool type_match(GeometricObject* pObj)	{return dynamic_cast<ConstrainedEdge*>(pObj) != NULL;}

		ConstrainedEdge() : m_pConstrainingObject(NULL), m_parentBaseObjectId(-1)	{}
		ConstrainedEdge(VertexBase* v1, VertexBase* v2) :
			m_pConstrainingObject(NULL),
			m_parentBaseObjectId(-1)
			{
				m_vertices[0] = v1;
				m_vertices[1] = v2;
			}

		ConstrainedEdge(const EdgeDescriptor& descriptor) :
			m_pConstrainingObject(NULL),
			m_parentBaseObjectId(-1)
			{
				m_vertices[0] = descriptor.vertex(0);
				m_vertices[1] = descriptor.vertex(1);
			}

		virtual ~ConstrainedEdge()	{}

		virtual GeometricObject* create_empty_instance() const	{return new ConstrainedEdge;}

		virtual int container_section() const	{return CSEDGE_CONSTRAINED_EDGE;}
		virtual ReferenceObjectID reference_object_id() const {return ROID_EDGE;}

		virtual bool is_constrained() const			{return true;}

	///	virtual refine. Returns pointers to EdgeBase.
	/**
	 * create 2 new constrained edges, connecting the original edges end-points with vrtNew.
	 * please note that the caller is responsible to resolve any conflicts arising with
	 * existing links of constrained/constraining objects.
	 * \sa EdgeBase::refine.
	 */
		virtual bool refine(std::vector<EdgeBase*>& vNewEdgesOut,
							VertexBase* newVertex,
							VertexBase** pSubstituteVrts = NULL);

//TODO:	Think about this method. It is not safe!
	///	non virtual refine. Returns pointers to ConstrainedEdge.
	/**
	 * create 2 new constrained edges, connecting the original edges end-points with vrtNew.
	 * please note that the caller is responsible to resolve any conflicts arising with
	 * existing links of constrained/constraining objects.
	 * \sa EdgeBase::refine
	 */
		bool refine(std::vector<ConstrainedEdge*>& vNewEdgesOut,
					VertexBase* newVertex,
					VertexBase** pSubstituteVrts = NULL);

		inline void set_constraining_object(GeometricObject* pObj)
		{
			m_pConstrainingObject = pObj;
			if(pObj)
				m_parentBaseObjectId = pObj->base_object_id();
		}

		inline GeometricObject* get_constraining_object()	{return m_pConstrainingObject;}

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
		GeometricObject*	m_pConstrainingObject;
		int					m_parentBaseObjectId;
};

template <>
class geometry_traits<ConstrainedEdge>
{
	public:
		typedef GenericGeometricObjectIterator<ConstrainedEdge*, EdgeBaseIterator>		iterator;
		typedef ConstGenericGeometricObjectIterator<ConstrainedEdge*, EdgeBaseIterator,
																ConstEdgeBaseIterator>	const_iterator;

		typedef EdgeDescriptor	Descriptor;
		typedef EdgeBase		geometric_base_object;

		enum
		{
			CONTAINER_SECTION = CSEDGE_CONSTRAINED_EDGE,
			BASE_OBJECT_ID = EDGE
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_EDGE;
};

typedef geometry_traits<ConstrainedEdge>::iterator 			ConstrainedEdgeIterator;
typedef geometry_traits<ConstrainedEdge>::const_iterator 	ConstConstrainedEdgeIterator;



////////////////////////////////////////////////////////////////////////
//	ConstrainingEdge
///	contains elements of type \sa HangingVertex and \sa ConstrainedEdge
/**
 * Edges of this type appear during hanging-vertex-refinement.
 * Treat with care.
 *
 * \ingroup lib_grid_geometric_objects
 */
class UG_API ConstrainingEdge : public EdgeBase
{
	friend class Grid;
	public:
		inline static bool type_match(GeometricObject* pObj)	{return dynamic_cast<ConstrainingEdge*>(pObj) != NULL;}

		ConstrainingEdge()	{}
		ConstrainingEdge(VertexBase* v1, VertexBase* v2)
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

		virtual ~ConstrainingEdge()	{}

		virtual GeometricObject* create_empty_instance() const	{return new ConstrainingEdge;}

		virtual int container_section() const	{return CSEDGE_CONSTRAINING_EDGE;}
		virtual ReferenceObjectID reference_object_id() const {return ROID_EDGE;}

		virtual bool is_constraining() const	{return true;}

	///	virtual refine. Returns pointers to EdgeBase.
	/**
	 * create 2 new constraining edges, connecting the original edges end-points with vrtNew.
	 * The refine has no effect on constrained objects. The user has to manually copy them.
	 * \sa EdgeBase::refine.
	 */
		virtual bool refine(std::vector<EdgeBase*>& vNewEdgesOut,
							VertexBase* newVertex,
							VertexBase** pSubstituteVrts = NULL);

//TODO:	Think about this method. It is not safe!
	///	non virtual refine. Returns pointers to ConstrainingEdge.
	/**
	 * create 2 new constraining edges, connecting the original edges end-points with vrtNew.
	 * The refine has no effect on constrained objects. The user has to manually copy them.
	 * \sa EdgeBase::refine.
	 */
		bool refine(std::vector<ConstrainingEdge*>& vNewEdgesOut,
						VertexBase* newVertex,
						VertexBase** pSubstituteVrts = NULL);


		inline void add_constrained_object(VertexBase* pObj)
			{
				UG_ASSERT(!is_constrained_object(pObj), "vertex is already constrained by this edge.");
					m_constrainedVertices.push_back(pObj);
			}

		inline void add_constrained_object(EdgeBase* pObj)
			{
				UG_ASSERT(!is_constrained_object(pObj), "edge is already constrained by this edge.");
					m_constrainedEdges.push_back(pObj);
			}

		inline bool is_constrained_object(VertexBase* vrt)
			{
				std::vector<VertexBase*>::iterator iter = find(m_constrainedVertices.begin(),
																m_constrainedVertices.end(), vrt);
				return iter != m_constrainedVertices.end();
			}

		inline bool is_constrained_object(EdgeBase* edge)
			{
				std::vector<EdgeBase*>::iterator iter = find(m_constrainedEdges.begin(),
																m_constrainedEdges.end(), edge);
				return iter != m_constrainedEdges.end();
			}

		inline void unconstrain_vertex(const VertexBase* vrt)
			{
				std::vector<VertexBase*>::iterator iter = find(m_constrainedVertices.begin(),
																m_constrainedVertices.end(), vrt);
				if(iter != m_constrainedVertices.end())
					m_constrainedVertices.erase(iter);
			}

		inline void unconstrain_edge(const EdgeBase* edge)
			{
				std::vector<EdgeBase*>::iterator iter = find(m_constrainedEdges.begin(),
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

		inline size_t num_constrained_vertices()	{return m_constrainedVertices.size();}
		inline size_t num_constrained_edges()		{return m_constrainedEdges.size();}

		VertexBase* constrained_vertex(size_t ind)
		{
			UG_ASSERT(ind < m_constrainedVertices.size(), "bad index");
			return m_constrainedVertices[ind];
		}

		EdgeBase* constrained_edge(size_t ind)
		{
			UG_ASSERT(ind < m_constrainedEdges.size(), "bad index");
			return m_constrainedEdges[ind];
		}

	protected:
		std::vector<VertexBase*>	m_constrainedVertices;
		std::vector<EdgeBase*>		m_constrainedEdges;
};

template <>
class geometry_traits<ConstrainingEdge>
{
	public:
		typedef GenericGeometricObjectIterator<ConstrainingEdge*, EdgeBaseIterator>			iterator;
		typedef ConstGenericGeometricObjectIterator<ConstrainingEdge*, EdgeBaseIterator,
																	ConstEdgeBaseIterator>	const_iterator;

		typedef EdgeDescriptor	Descriptor;
		typedef EdgeBase		geometric_base_object;

		enum
		{
			CONTAINER_SECTION = CSEDGE_CONSTRAINING_EDGE,
			BASE_OBJECT_ID = EDGE
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_EDGE;
};

typedef geometry_traits<ConstrainingEdge>::iterator 		ConstrainingEdgeIterator;
typedef geometry_traits<ConstrainingEdge>::const_iterator 	ConstConstrainingEdgeIterator;

}//	end of namespace

#endif
