// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 23.12.2011 (m,d,y)

#ifndef __H__UG__geometric_objects_2d__
#define __H__UG__geometric_objects_2d__

#include "../grid/grid.h"
#include "common/math/ugmath.h"
#include "common/assert.h"
#include "geometric_objects_0d.h"
#include "geometric_objects_1d.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
///	These numbers define where in the face-section-container a face will be stored.
/**	The order of the constants must not be changed! Algorithms may exist that rely on it.*/
enum FaceContainerSections
{
	CSFACE_NONE = -1,
	CSFACE_TRIANGLE = 0,
	CSFACE_QUADRILATERAL = 1,
	CSFACE_CONSTRAINED_TRIANGLE = 2,
	CSFACE_CONSTRAINED_QUADRILATERAL = 3,
	CSFACE_CONSTRAINING_TRIANGLE = 4,
	CSFACE_CONSTRAINING_QUADRILATERAL = 5,

	CSFACE_USER	// always last
};


////////////////////////////////////////////////////////////////////////////////
//	NORMAL FACES
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//	TriangleDescriptor
///	only used to initialize a triangle. for all other tasks you should use FaceDescriptor.
class UG_API TriangleDescriptor
{
	public:
		TriangleDescriptor()	{}
		TriangleDescriptor(const TriangleDescriptor& td);
		TriangleDescriptor(VertexBase* v1, VertexBase* v2, VertexBase* v3);

		inline uint num_vertices() const					{return 3;}
		inline void set_vertex(uint index, VertexBase* v)	{m_vertex[index] = v;}
		inline VertexBase* vertex(uint index) const			{return m_vertex[index];}

	protected:
		VertexBase*	m_vertex[3];
};



////////////////////////////////////////////////////////////////////////
//	CustomTriangle
///	Concrete types share this base-type. It is not intended for direct use.
/**
 * BaseClass has to be derived from Face (or simply should be Face).
 * The ConcreteTriangleType is used in methods like refine, etc. as the type
 * of newly created objects.
 */
template <class ConcreteTriangleType, class BaseClass>
class UG_API CustomTriangle : public BaseClass
{
	public:
		CustomTriangle()	{}
		CustomTriangle(const TriangleDescriptor& td);
		CustomTriangle(VertexBase* v1, VertexBase* v2, VertexBase* v3);

		virtual GeometricObject* create_empty_instance() const	{return new ConcreteTriangleType;}
		virtual ReferenceObjectID reference_object_id() const {return ROID_TRIANGLE;}

		virtual VertexBase* vertex(uint index) const	{return m_vertices[index];}
		virtual Face::ConstVertexArray vertices() const		{return m_vertices;}
		virtual size_t num_vertices() const	{return 3;}

		virtual EdgeDescriptor edge_desc(int index) const
			{return EdgeDescriptor(m_vertices[index], m_vertices[(index+1) % 3]);}

		virtual void edge_desc(int index, EdgeDescriptor& edOut) const
			{edOut.set_vertices(m_vertices[index], m_vertices[(index+1) % 3]);}

		virtual std::pair<GeometricBaseObject, int> get_opposing_object(VertexBase* vrt) const;

	///	Refines a Triangle by inserting new vertices. \sa Face::refine.
		virtual bool refine(std::vector<Face*>& vNewFacesOut,
							VertexBase** newFaceVertexOut,
							VertexBase** newEdgeVertices,
							VertexBase* newFaceVertex = NULL,
							VertexBase** pSubstituteVertices = NULL);

		virtual bool collapse_edge(std::vector<Face*>& vNewFacesOut,
								int edgeIndex, VertexBase* newVertex,
								VertexBase** pSubstituteVertices = NULL);

		virtual bool collapse_edges(std::vector<Face*>& vNewFacesOut,
								std::vector<VertexBase*>& vNewEdgeVertices,
								VertexBase** pSubstituteVertices = NULL);

//	BEGIN Depreciated
		virtual void create_faces_by_edge_split(int splitEdgeIndex,
							VertexBase* newVertex,
							std::vector<Face*>& vNewFacesOut,
							VertexBase** pSubstituteVertices = NULL);

	protected:
		virtual void set_vertex(uint index, VertexBase* pVrt)	{m_vertices[index] = pVrt;}

	protected:
		VertexBase* m_vertices[3];
};



////////////////////////////////////////////////////////////////////////
//	Triangle
///	the most simple form of a face
/**
 *
 * \ingroup lib_grid_geometric_objects
 */
class UG_API Triangle : public CustomTriangle<Triangle, Face>
{
	typedef CustomTriangle<Triangle, Face> BaseClass;
	public:
		inline static bool type_match(GeometricObject* pObj)	{return dynamic_cast<Triangle*>(pObj) != NULL;}

		Triangle() : BaseClass()	{}
		Triangle(const TriangleDescriptor& td) : BaseClass(td)	{}
		Triangle(VertexBase* v1, VertexBase* v2, VertexBase* v3) : BaseClass(v1, v2, v3)	{}

		virtual int container_section() const	{return CSFACE_TRIANGLE;}

	protected:
		virtual EdgeBase* create_edge(int index)
			{
				return new Edge(m_vertices[index], m_vertices[(index+1) % 3]);
			}
};

template <>
class geometry_traits<Triangle>
{
	public:
		typedef GenericGeometricObjectIterator<Triangle*, FaceIterator>			iterator;
		typedef ConstGenericGeometricObjectIterator<Triangle*, FaceIterator,
															ConstFaceIterator>	const_iterator;

		typedef TriangleDescriptor Descriptor;	///< Faces can't be created directly
		typedef Face	geometric_base_object;

		enum
		{
			CONTAINER_SECTION = CSFACE_TRIANGLE,
			BASE_OBJECT_ID = FACE
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_TRIANGLE;
};

typedef geometry_traits<Triangle>::iterator			TriangleIterator;
typedef geometry_traits<Triangle>::const_iterator	ConstTriangleIterator;



////////////////////////////////////////////////////////////////////////
//	QuadrilateralDescriptor
///	only used to initialize a quadrilateral. for all other tasks you should use FaceDescriptor.
class UG_API QuadrilateralDescriptor
{
	public:
		QuadrilateralDescriptor()	{}
		QuadrilateralDescriptor(const QuadrilateralDescriptor& qd);
		QuadrilateralDescriptor(VertexBase* v1, VertexBase* v2, VertexBase* v3, VertexBase* v4);

		inline uint num_vertices() const					{return 4;}
		inline void set_vertex(uint index, VertexBase* v)	{m_vertex[index] = v;}
		inline VertexBase* vertex(uint index) const			{return m_vertex[index];}

	protected:
		VertexBase*	m_vertex[4];
};



////////////////////////////////////////////////////////////////////////
//	CustomQuadrilateral
///	Concrete types share this base-type. It is not intended for direct use.
/**
 * BaseClass has to be derived from Face (or simply should be Face).
 * The ConcreteQuadrilateralType is used in methods like refine, etc. as the type
 * of newly created objects.
 */
template <class ConcreteQuadrilateralType, class BaseClass>
class UG_API CustomQuadrilateral : public BaseClass
{
	public:
		using Face::ConstVertexArray;

		CustomQuadrilateral()	{}
		CustomQuadrilateral(const QuadrilateralDescriptor& qd);
		CustomQuadrilateral(VertexBase* v1, VertexBase* v2,
							VertexBase* v3, VertexBase* v4);

		virtual GeometricObject* create_empty_instance() const	{return new ConcreteQuadrilateralType;}
		virtual ReferenceObjectID reference_object_id() const {return ROID_QUADRILATERAL;}

		virtual VertexBase* vertex(uint index) const	{return m_vertices[index];}
		virtual Face::ConstVertexArray vertices() const		{return m_vertices;}
		virtual size_t num_vertices() const	{return 4;}

		virtual EdgeDescriptor edge_desc(int index) const
			{return EdgeDescriptor(m_vertices[index], m_vertices[(index+1) % 4]);}

		virtual void edge_desc(int index, EdgeDescriptor& edOut) const
			{edOut.set_vertices(m_vertices[index], m_vertices[(index+1) % 4]);}


	///	fills the edge-descriptor with the edge that lies opposed to the specified one
	/**	If the specified edge is not part of the face, false is returned.*/
		virtual bool get_opposing_side(EdgeVertices* e, EdgeDescriptor& edOut) const;

		virtual std::pair<GeometricBaseObject, int> get_opposing_object(VertexBase* vrt) const;

	///	Refines a Quadrilateral by inserting new vertices. \sa Face::refine.
		virtual bool refine(std::vector<Face*>& vNewFacesOut,
							VertexBase** newFaceVertexOut,
							VertexBase** newEdgeVertices,
							VertexBase* newFaceVertex = NULL,
							VertexBase** pSubstituteVertices = NULL);

		virtual bool collapse_edge(std::vector<Face*>& vNewFacesOut,
								int edgeIndex, VertexBase* newVertex,
								VertexBase** pSubstituteVertices = NULL);

		virtual bool collapse_edges(std::vector<Face*>& vNewFacesOut,
								std::vector<VertexBase*>& vNewEdgeVertices,
								VertexBase** pSubstituteVertices = NULL);

//	BEGIN Depreciated
		virtual void create_faces_by_edge_split(int splitEdgeIndex,
							VertexBase* newVertex,
							std::vector<Face*>& vNewFacesOut,
							VertexBase** pSubstituteVertices = NULL);

	protected:
		virtual void set_vertex(uint index, VertexBase* pVrt)	{m_vertices[index] = pVrt;}

	protected:
		VertexBase* m_vertices[4];
};



////////////////////////////////////////////////////////////////////////
//	Quadrilateral
///	a face with four points.
/**
 * \ingroup lib_grid_geometric_objects
 */
class UG_API Quadrilateral : public CustomQuadrilateral<Quadrilateral, Face>
{
	public:
		typedef CustomQuadrilateral<Quadrilateral, Face> BaseClass;

		inline static bool type_match(GeometricObject* pObj)	{return dynamic_cast<Quadrilateral*>(pObj) != NULL;}

		Quadrilateral()	{}
		Quadrilateral(const QuadrilateralDescriptor& td) : BaseClass(td)	{}
		Quadrilateral(VertexBase* v1, VertexBase* v2,
					  VertexBase* v3, VertexBase* v4) : BaseClass(v1, v2, v3, v4)	{}

		virtual int container_section() const	{return CSFACE_QUADRILATERAL;}

	protected:
		virtual EdgeBase* create_edge(int index)
		{
			return new Edge(m_vertices[index], m_vertices[(index+1) % 4]);
		}
};

template <>
class geometry_traits<Quadrilateral>
{
	public:
		typedef GenericGeometricObjectIterator<Quadrilateral*, FaceIterator>		iterator;
		typedef ConstGenericGeometricObjectIterator<Quadrilateral*, FaceIterator,
																ConstFaceIterator>	const_iterator;

		typedef QuadrilateralDescriptor Descriptor;	///< Faces can't be created directly
		typedef Face	geometric_base_object;

		enum
		{
			CONTAINER_SECTION = CSFACE_QUADRILATERAL,
			BASE_OBJECT_ID = FACE
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_QUADRILATERAL;
};

typedef geometry_traits<Quadrilateral>::iterator		QuadrilateralIterator;
typedef geometry_traits<Quadrilateral>::const_iterator	ConstQuadrilateralIterator;



////////////////////////////////////////////////////////////////////////////////
//	CONSTRAINED FACES
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//	ConstrainedFace
///	This class stores the constraining object.
/**
 * Please note, that the user is has to link and unlink constraining
 * objects manually.
 */
class UG_API ConstrainedFace : public Face
{
	public:
		inline static bool type_match(GeometricObject* pObj)	{return dynamic_cast<ConstrainedFace*>(pObj) != NULL;}

		ConstrainedFace() : m_pConstrainingObject(NULL)			{}
		virtual ~ConstrainedFace()	{}

		inline void set_constraining_object(GeometricObject* pObj)	{m_pConstrainingObject = pObj;}
		inline GeometricObject* get_constraining_object()			{return m_pConstrainingObject;}

		virtual bool is_constrained() const							{return true;}

	protected:
		GeometricObject*	m_pConstrainingObject;
};


////////////////////////////////////////////////////////////////////////
//	ConstrainedTriangle
///	a triangle constrained by another object.
/**
 * \ingroup lib_grid_geometric_objects
 */
class UG_API ConstrainedTriangle : public CustomTriangle<ConstrainedTriangle, ConstrainedFace>
{
	typedef CustomTriangle<ConstrainedTriangle, ConstrainedFace> BaseTriangle;

	public:
		inline static bool type_match(GeometricObject* pObj)	{return dynamic_cast<ConstrainedTriangle*>(pObj) != NULL;}

		ConstrainedTriangle() :
			BaseTriangle()	{}

		ConstrainedTriangle(const TriangleDescriptor& td) :
			BaseTriangle(td)	{}

		ConstrainedTriangle(VertexBase* v1, VertexBase* v2, VertexBase* v3) :
			BaseTriangle(v1, v2, v3)	{}

		virtual int container_section() const	{return CSFACE_CONSTRAINED_TRIANGLE;}

	protected:
		virtual EdgeBase* create_edge(int index)
			{
				return new ConstrainedEdge(m_vertices[index], m_vertices[(index+1) % 3]);
			}
};

template <>
class geometry_traits<ConstrainedTriangle>
{
	public:
		typedef GenericGeometricObjectIterator<ConstrainedTriangle*, FaceIterator>		iterator;
		typedef ConstGenericGeometricObjectIterator<ConstrainedTriangle*, FaceIterator,
																	ConstFaceIterator>	const_iterator;

		typedef TriangleDescriptor Descriptor;	///< Faces can't be created directly
		typedef Face	geometric_base_object;

		enum
		{
			CONTAINER_SECTION = CSFACE_CONSTRAINED_TRIANGLE,
			BASE_OBJECT_ID = FACE
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_TRIANGLE;
};

typedef geometry_traits<ConstrainedTriangle>::iterator			ConstrainedTriangleIterator;
typedef geometry_traits<ConstrainedTriangle>::const_iterator	ConstConstrainedTriangleIterator;



////////////////////////////////////////////////////////////////////////
//	ConstrainedQuadrilateral
///	a quadrilateral constrained by another object.
/**
 * \ingroup lib_grid_geometric_objects
 */
class UG_API ConstrainedQuadrilateral : public CustomQuadrilateral<ConstrainedQuadrilateral, ConstrainedFace>
{
	typedef CustomQuadrilateral<ConstrainedQuadrilateral, ConstrainedFace> BaseClass;

	public:
		inline static bool type_match(GeometricObject* pObj)	{return dynamic_cast<ConstrainedQuadrilateral*>(pObj) != NULL;}

		ConstrainedQuadrilateral() : BaseClass()	{}
		ConstrainedQuadrilateral(const QuadrilateralDescriptor& qd) : BaseClass(qd)	{}
		ConstrainedQuadrilateral(VertexBase* v1, VertexBase* v2,
								 VertexBase* v3, VertexBase* v4) : BaseClass(v1, v2, v3, v4)	{}

		virtual int container_section() const	{return CSFACE_CONSTRAINED_QUADRILATERAL;}

	protected:
		virtual EdgeBase* create_edge(int index)
			{
				return new ConstrainedEdge(m_vertices[index], m_vertices[(index+1) % 4]);
			}
};


template <>
class geometry_traits<ConstrainedQuadrilateral>
{
	public:
		typedef GenericGeometricObjectIterator<ConstrainedQuadrilateral*, FaceIterator>	iterator;
		typedef ConstGenericGeometricObjectIterator<ConstrainedQuadrilateral*,
													FaceIterator, ConstFaceIterator>	const_iterator;

		typedef QuadrilateralDescriptor Descriptor;	///< Faces can't be created directly
		typedef Face	geometric_base_object;

		enum
		{
			CONTAINER_SECTION = CSFACE_CONSTRAINED_QUADRILATERAL,
			BASE_OBJECT_ID = FACE
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_QUADRILATERAL;
};

typedef geometry_traits<ConstrainedQuadrilateral>::iterator			ConstrainedQuadrilateralIterator;
typedef geometry_traits<ConstrainedQuadrilateral>::const_iterator	ConstConstrainedQuadrilateralIterator;



////////////////////////////////////////////////////////////////////////////////
//	CONSTRAINING FACES
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//	ConstrainingFace
///	This class is used to store constrained geometric objects.
/**
 * Please note, that the user is has to link and unlink constrained
 * objects manually.
 */
class UG_API ConstrainingFace : public Face
{

	public:
		inline static bool type_match(GeometricObject* pObj)	{return dynamic_cast<ConstrainingFace*>(pObj) != NULL;}

		virtual ~ConstrainingFace()	{}

		virtual bool is_constraining() const					{return true;}

		inline void add_constrained_object(VertexBase* pObj)
			{
				UG_ASSERT(!is_constrained_object(pObj), "vertex is already registered at constraining face");
					m_constrainedVertices.push_back(pObj);
			}

		inline void add_constrained_object(EdgeBase* pObj)
			{
				UG_ASSERT(!is_constrained_object(pObj), "edge is already registered at constraining face");
					m_constrainedEdges.push_back(pObj);
			}

		inline void add_constrained_object(Face* pObj)
			{
				UG_ASSERT(!is_constrained_object(pObj), "face is already registered at constraining face");
					m_constrainedFaces.push_back(pObj);
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

		inline bool is_constrained_object(Face* face)
			{
				std::vector<Face*>::iterator iter = find(m_constrainedFaces.begin(),
														m_constrainedFaces.end(), face);
				return iter != m_constrainedFaces.end();
			}

		inline void unconstrain_object(const VertexBase* vrt)
			{
				std::vector<VertexBase*>::iterator iter = find(m_constrainedVertices.begin(),
															 m_constrainedVertices.end(), vrt);
				if(iter != m_constrainedVertices.end())
					m_constrainedVertices.erase(iter);
			}

		inline void unconstrain_object(const EdgeBase* edge)
			{
				std::vector<EdgeBase*>::iterator iter = find(m_constrainedEdges.begin(),
															m_constrainedEdges.end(), edge);
				if(iter != m_constrainedEdges.end())
					m_constrainedEdges.erase(iter);
			}

		inline void unconstrain_object(const Face* face)
			{
				std::vector<Face*>::iterator iter = find(m_constrainedFaces.begin(),
														m_constrainedFaces.end(), face);
				if(iter != m_constrainedFaces.end())
					m_constrainedFaces.erase(iter);
			}

		inline void clear_constrained_vertices()	{m_constrainedVertices.clear();}
		inline void clear_constrained_edges()		{m_constrainedEdges.clear();}
		inline void clear_constrained_faces()		{m_constrainedFaces.clear();}
		inline void clear_constrained_objects()
			{
				clear_constrained_vertices();
				clear_constrained_edges();
				clear_constrained_faces();
			}

		inline size_t num_constrained_vertices()	{return m_constrainedVertices.size();}
		inline size_t num_constrained_edges()		{return m_constrainedEdges.size();}
		inline size_t num_constrained_faces()		{return m_constrainedFaces.size();}

		inline VertexBase* constrained_vertex(size_t ind)
			{
				UG_ASSERT(ind < m_constrainedVertices.size(), "bad index.");
				return m_constrainedVertices[ind];
			}

		inline EdgeBase* constrained_edge(size_t ind)
			{
				UG_ASSERT(ind < m_constrainedEdges.size(), "bad index.");
				return m_constrainedEdges[ind];
			}

		inline Face* constrained_face(size_t ind)
			{
				UG_ASSERT(ind < m_constrainedFaces.size(), "bad index.");
				return m_constrainedFaces[ind];
			}

	protected:
		std::vector<VertexBase*>	m_constrainedVertices;
		std::vector<EdgeBase*>		m_constrainedEdges;
		std::vector<Face*>			m_constrainedFaces;
};



////////////////////////////////////////////////////////////////////////
//	ConstrainingTriangle
///	a triangle constraining other objects.
/**
 * \ingroup lib_grid_geometric_objects
 */
class UG_API ConstrainingTriangle : public CustomTriangle<ConstrainingTriangle, ConstrainingFace>
{
	typedef CustomTriangle<ConstrainingTriangle, ConstrainingFace> BaseTriangle;

	public:
		inline static bool type_match(GeometricObject* pObj)	{return dynamic_cast<ConstrainingTriangle*>(pObj) != NULL;}

		ConstrainingTriangle() :
			BaseTriangle()	{reserve_memory();}
		ConstrainingTriangle(const TriangleDescriptor& td) :
			BaseTriangle(td)	{reserve_memory();}
		ConstrainingTriangle(VertexBase* v1, VertexBase* v2, VertexBase* v3) :
			BaseTriangle(v1, v2, v3)	{reserve_memory();}

		virtual int container_section() const	{return CSFACE_CONSTRAINING_TRIANGLE;}

	protected:
		void reserve_memory()
			{
				m_constrainedEdges.reserve(3);
				m_constrainedFaces.reserve(4);
			}

		virtual EdgeBase* create_edge(int index)
			{
				return new Edge(m_vertices[index], m_vertices[(index+1) % 3]);
			}
};

template <>
class geometry_traits<ConstrainingTriangle>
{
	public:
		typedef GenericGeometricObjectIterator<ConstrainingTriangle*, FaceIterator>		iterator;
		typedef ConstGenericGeometricObjectIterator<ConstrainingTriangle*, FaceIterator,
																	ConstFaceIterator>	const_iterator;

		typedef TriangleDescriptor Descriptor;	///< Faces can't be created directly
		typedef Face	geometric_base_object;

		enum
		{
			CONTAINER_SECTION = CSFACE_CONSTRAINING_TRIANGLE,
			BASE_OBJECT_ID = FACE
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_TRIANGLE;
};

typedef geometry_traits<ConstrainingTriangle>::iterator			ConstrainingTriangleIterator;
typedef geometry_traits<ConstrainingTriangle>::const_iterator	ConstConstrainingTriangleIterator;



////////////////////////////////////////////////////////////////////////
//	ConstrainingQuadrilateral
///	a quadrilateral constraining other objects.
/**
 * \ingroup lib_grid_geometric_objects
 */
class UG_API ConstrainingQuadrilateral : public CustomQuadrilateral<ConstrainingQuadrilateral, ConstrainingFace>
{
	typedef CustomQuadrilateral<ConstrainingQuadrilateral, ConstrainingFace> BaseClass;

	public:
		inline static bool type_match(GeometricObject* pObj)	{return dynamic_cast<ConstrainingQuadrilateral*>(pObj) != NULL;}

		ConstrainingQuadrilateral() :
			BaseClass()	{reserve_memory();}
		ConstrainingQuadrilateral(const QuadrilateralDescriptor& qd) :
			BaseClass(qd)	{reserve_memory();}
		ConstrainingQuadrilateral(VertexBase* v1, VertexBase* v2,
								  VertexBase* v3, VertexBase* v4) :
			BaseClass(v1, v2, v3, v4)	{reserve_memory();}

		virtual int container_section() const	{return CSFACE_CONSTRAINING_QUADRILATERAL;}

	protected:
		void reserve_memory()
			{
				m_constrainedVertices.reserve(1);
				m_constrainedEdges.reserve(4);
				m_constrainedFaces.reserve(4);
			}

		virtual EdgeBase* create_edge(int index)
			{
				return new Edge(m_vertices[index], m_vertices[(index+1) % 4]);
			}
};

template <>
class geometry_traits<ConstrainingQuadrilateral>
{
	public:
		typedef GenericGeometricObjectIterator<ConstrainingQuadrilateral*, FaceIterator>	iterator;
		typedef ConstGenericGeometricObjectIterator<ConstrainingQuadrilateral*,
														FaceIterator, ConstFaceIterator>	const_iterator;

		typedef QuadrilateralDescriptor Descriptor;	///< Faces can't be created directly
		typedef Face	geometric_base_object;

		enum
		{
			CONTAINER_SECTION = CSFACE_CONSTRAINING_QUADRILATERAL,
			BASE_OBJECT_ID = FACE
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_QUADRILATERAL;
};

typedef geometry_traits<ConstrainingQuadrilateral>::iterator		ConstrainingQuadrilateralIterator;
typedef geometry_traits<ConstrainingQuadrilateral>::const_iterator	ConstConstrainingQuadrilateralIterator;


}//	end of namespace

#endif
