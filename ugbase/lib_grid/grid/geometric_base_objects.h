//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d09

#ifndef __H__LIB_GRID__GEOMETRIC_BASE_OBJECTS__
#define __H__LIB_GRID__GEOMETRIC_BASE_OBJECTS__

#include <list>
#include <cassert>
#include "common/types.h"
#include "common/util/attachment_pipe.h"
#include "common/util/hash.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////////////////////////////
//	Predeclaration of the Grid-type.
class Grid;

////////////////////////////////////////////////////////////////////////
//	GeometricBaseObject
///	enumeration of the GeometricBaseObjects that make up a grid.
enum GeometricBaseObject
{
	VERTEX = 0,
	EDGE,
	FACE,
	VOLUME,
	NUM_GEOMETRIC_BASE_OBJECTS		//always last!!!
};

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	Reference-Object IDs
///	these ids are used to identify the shape of a geometric object.
enum ReferenceObjectID
{
	ROID_INVALID = -1,
	ROID_VERTEX,
	ROID_EDGE,
	ROID_TRIANGLE,
	ROID_QUADRILATERAL,
	ROID_TETRAHEDRON,
	ROID_HEXAHEDRON,
	ROID_PRISM,
	ROID_PYRAMID,
	NUM_REFERENCE_OBJECTS
};


////////////////////////////////////////////////////////////////////////
//	Predeclaration of geometric objects.
class GeometricObject;	//	geometric base object
class VertexBase;		//	base for all 0-dimensional grid objects.
class EdgeBase;			//	base for all 1-dimensional grid objects.
class Face;				//	base for all 2-dimensional grid objects.
class Volume;			//	base for all 3-dimensional grid objects.

class EdgeDescriptor;	//	describes an edge.
class FaceDescriptor;	//	describes a face.
class VolumeDescriptor;	//	describes a volume.

class EdgeVertices;		//	manages the vertices of an edge. Base for EdgeBase and EdgeDescriptor.
class FaceVertices;		//	manages the vertices of a face. Base for Face and FaceDescriptor.
class VolumeVertices;	//	manages the vertices of a volume. Base for Volume and VolumeDescriptor.

//	pointer-types. Primarily required for template-specializations.
typedef VertexBase*	PVertexBase;
typedef EdgeBase*		PEdgeBase;
typedef Face*			PFace;
typedef Volume*		PVolume;

typedef EdgeDescriptor*	PEdgeDescriptor;
typedef FaceDescriptor*	PFaceDescriptor;
typedef VolumeDescriptor*	PVolumeDescriptor;

typedef EdgeVertices*		PEdgeVertices;
typedef FaceVertices*		PFaceVertices;
typedef VolumeVertices*	PVolumeVertices;



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	GeometricObjectContainer
///	Declaration of the container that will hold all geometric objects.
typedef std::list<GeometricObject*> 		GeometricObjectContainer;

////////////////////////////////////////////////////////////////////////
//	GeometricObjectIterator
///	This Iterator will be used as base-class for iterators of specialized geometric objects.
typedef GeometricObjectContainer::iterator			GeometricObjectIterator;
typedef GeometricObjectContainer::const_iterator	ConstGeometricObjectIterator;

////////////////////////////////////////////////////////////////////////
//	The geometry_traits. This class can be specialized by each element-type.
/**
 * In order to use a custom geometric object with libGrid, you have to
 * supply a specialization for geometry_traits.
 * Specializations have to specify the following types and methods:
 *
 * MANDATORY:
 * Types:
 * - geometric_base_object:		the geometric object from which TElem derives.
 * 					has to be either VertexBase, EdgeBase, Face or Volume.
 * - iterator:		An iterator that iterates over ElementContainer<BaseClass>
 * 					and which has a constructor that takes
 * 					ElementContainer<BaseClass>::iterator as an argument.
 * 					casts should be checked when in DEBUG mode!
 *
 * constants:
 * - SHARED_PIPE_SECTION: This constant should hold the pipes section in which your objects should be placed after creation, starting from 0. See the Grid-documentation for more information.
 * - BASE_OBJECT_TYPE_ID: Has to hold one of the GeometricBaseObject constants, or -1.
 *
 * OPTIONAL:
 * Types:
 * - Descriptor:	a class which can be passed to the constructor of the element.
 */
template <class TElem>
class geometry_traits
{};

////////////////////////////////////////////////////////////////////////////////////////////////
//	GenericGeometricObjectIterator
///	Use this class as a tool to create iterators to your own geometric objects.
/**
 *
 */
template <class TValue, class TBaseIterator = GeometricObjectIterator>
class GenericGeometricObjectIterator : public TBaseIterator
{
	friend class Grid;
	template <class TIterDest, class TIterSrc> friend TIterDest iterator_cast(const TIterSrc& iter);

	public:
		typedef TValue	value_type;

	public:
		GenericGeometricObjectIterator()	{}

		GenericGeometricObjectIterator(const GenericGeometricObjectIterator& iter) :
			TBaseIterator(iter)	{}

		inline TValue& operator* ()	{return (TValue&)(TBaseIterator::operator*());}
		inline const TValue& operator* () const	{return (TValue&)(TBaseIterator::operator*());}

	protected:
		GenericGeometricObjectIterator(const GeometricObjectIterator& iter) :
			TBaseIterator(iter)	{}
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	ConstGenericGeometricObjectIterator
///	Use this class as a tool to create const_iterators to your own geometric objects.
/**
 *
 */
template <class TValue, class TBaseIterator = ConstGeometricObjectIterator>
class ConstGenericGeometricObjectIterator : public TBaseIterator
{
	friend class Grid;
	template <class TIterDest, class TIterSrc> friend TIterDest iterator_cast(const TIterSrc& iter);

	public:
		typedef TValue	value_type;

	public:
		ConstGenericGeometricObjectIterator()	{}

		ConstGenericGeometricObjectIterator(const ConstGenericGeometricObjectIterator& iter) :
			TBaseIterator(iter)	{}

		inline const TValue& operator* () const	{return (TValue&)(TBaseIterator::operator*());}

	protected:
		ConstGenericGeometricObjectIterator(const GeometricObjectIterator& iter) :
			TBaseIterator(iter)	{}

		ConstGenericGeometricObjectIterator(const ConstGeometricObjectIterator& iter) :
			TBaseIterator(iter)	{}
};


////////////////////////////////////////////////////////////////////////
//	iterator_cast
///	You should avoid casting whenever possible!
template <class TIterDest, class TIterSrc>
inline TIterDest
iterator_cast(const TIterSrc& iter)
{
	return TIterDest(iter);
}

/**
 * \brief Geometric objects are the building blocks of a grid.
 *
 * \defgroup lib_grid_geometric_objects geometric objects
 * \ingroup lib_grid
 */
////////////////////////////////////////////////////////////////////////
//	GeometricObject
///	The base class for all geometric objects, such as vertices, edges, faces, volumes, ...
/**
 * In order to be used by libGrid, all derivatives of GeometricObject
 * have to specialize geometry_traits<GeomObjectType>.
 *
 * \ingroup lib_grid_geometric_objects
 */
class GeometricObject
{
	friend class Grid;
	friend class attachment_traits<GeometricObject*, Grid>;
	public:
		virtual ~GeometricObject()	{}

	///	create an instance of the derived type
	/**	Make sure to overload this method in derivates of this class!*/
		virtual GeometricObject* create_empty_instance() const {return NULL;}

		virtual int shared_pipe_section() const = 0;
		virtual int base_object_type_id() const = 0;//	This method probably shouldn't be there!
	/**
	 * A reference object represents a class of geometric objects.
	 * Tetrahedrons, Triangles etc are such classes.
	 * Reference ids should be defined in the file in which concrete geometric objects are defined.
	 */
		virtual int reference_object_id() const = 0;///	returns the id of the reference-object.

	protected:
		uint						m_gridDataIndex;//	index to grid-attached data.
		GeometricObjectIterator		m_entryIter;//	used for fast list-access.
};

template <>
class geometry_traits<GeometricObject>
{
	public:
		typedef GeometricObjectIterator			iterator;
		typedef ConstGeometricObjectIterator	const_iterator;

		enum
		{
			SHARED_PIPE_SECTION = -1,
			BASE_OBJECT_TYPE_ID = -1
		};
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	VertexBase
///	Base-class for all vertex-types
/**
 * Vertices are required in any grid.
 * They are the geometric objects of lowest dimension.
 * All other geometric objects of higher dimension reference vertices.
 *
 * \ingroup lib_grid_geometric_objects
 */
class VertexBase : public GeometricObject
{
	friend class Grid;
	public:
		inline static bool type_match(GeometricObject* pObj)	{return dynamic_cast<VertexBase*>(pObj) != NULL;}

		virtual ~VertexBase()	{}

		virtual int shared_pipe_section() const	{return -1;}
		virtual int base_object_type_id() const	{return VERTEX;}
		virtual int reference_object_id() const	{return -1;}

	///	returns a value that can be used for hashing.
	/**	this value is unique per grid.
	 *	It is only set correctly if the vertex has been created
	 *	by the grid or has been properly registered at it.*/
		inline uint32 get_hash_value() const	{return m_hashValue;}

	protected:
		uint32	m_hashValue;//	a unique value for each vertex in a grid.
};

template <>
class geometry_traits<VertexBase>
{
	public:
		typedef GenericGeometricObjectIterator<VertexBase*>			iterator;
		typedef ConstGenericGeometricObjectIterator<VertexBase*>	const_iterator;

		typedef VertexBase	geometric_base_object;

		enum
		{
			SHARED_PIPE_SECTION = -1,
			BASE_OBJECT_TYPE_ID = VERTEX
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_VERTEX;
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	EdgeVertices
///	holds the vertices of an EdgeBase or an EdgeDescriptor.
/**	Please note that this class does not have a virtual destructor.*/
class EdgeVertices
{
	friend class Grid;
	public:
		inline VertexBase* vertex(uint index) const	{return m_vertices[index];}
		inline uint num_vertices() const			{return 2;}	// this method is supplied to allow the use of EdgeBase in template-methods that require a num_vertices() method.

	//	compatibility with std::vector for some template routines
	///	returns the number of vertices.
		inline size_t size() const	{return 2;}
	///	returns the i-th vertex.
		VertexBase* operator[](uint index) const {return m_vertices[index];}

	protected:
		inline void assign_edge_vertices(const EdgeVertices& ev)
		{
			m_vertices[0] = ev.m_vertices[0];
			m_vertices[1] = ev.m_vertices[1];
		}

	protected:
		VertexBase*	m_vertices[2];
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	EdgeBase
///	Base-class for edges
/**
 * EdgeBase is the base class of all 1-dimensional geometric objects.
 * Edges connect two vertices.
 *
 * \ingroup lib_grid_geometric_objects
 */
class EdgeBase : public GeometricObject, public EdgeVertices
{
	friend class Grid;
	public:
		// lower dimensional Base Object
		typedef VertexBase lower_dim_base_object;

	public:
		inline static bool type_match(GeometricObject* pObj)	{return dynamic_cast<EdgeBase*>(pObj) != NULL;}

		virtual ~EdgeBase()	{}

		virtual int shared_pipe_section() const	{return -1;}
		virtual int base_object_type_id() const	{return EDGE;}
		virtual int reference_object_id() const	{return -1;}

		virtual int num_sides() const {return 2;}
	/**
	 * create 2 new edges, connecting the original edges end-points with vrtNew.
	 * Newly created edges have to be registered at a grid manually by the caller.
	 * If the caller does not register the edges in vGeomOut at a grid, he is
	 * responsible to free the associated memory (delete each element in vNewEdgesOut).
	 * Please note that refining an edge using this method does not automatically
	 * refine associated elements.
	 * Be sure to store the new edges in the right order. vNewEdgesOut should contain
	 * the edge connecting vertex(0) and newVertex first.
	 *
	 * You may pass an array of 2 vertices to pSubstituteVrts. If you do so, Those
	 * vertices will be used instead of the original ones.
	 */
		virtual bool refine(std::vector<EdgeBase*>& vNewEdgesOut,
											VertexBase* newVertex,
											VertexBase** pSubstituteVrts = NULL)	{return false;}

	protected:
		inline void set_vertex(uint index, VertexBase* pVrt)	{m_vertices[index] = pVrt;}
};

template <>
class geometry_traits<EdgeBase>
{
	public:
		typedef GenericGeometricObjectIterator<EdgeBase*>		iterator;
		typedef ConstGenericGeometricObjectIterator<EdgeBase*>	const_iterator;

		typedef EdgeBase	geometric_base_object;

		enum
		{
			SHARED_PIPE_SECTION = -1,
			BASE_OBJECT_TYPE_ID = EDGE
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_EDGE;
};


////////////////////////////////////////////////////////////////////////////////////////////////
//	EdgeDescriptor
///	Can be used to store information about an edge and to construct an edge.
class EdgeDescriptor : public EdgeVertices
{
	public:
		EdgeDescriptor();
		EdgeDescriptor(const EdgeDescriptor& ed);
		EdgeDescriptor(VertexBase* vrt1, VertexBase* vrt2);

		EdgeDescriptor& operator = (const EdgeDescriptor& ed);

		inline void set_vertex(uint index, VertexBase* vrt)	{m_vertices[index] = vrt;}
		inline void set_vertices(VertexBase* vrt1, VertexBase* vrt2)
			{
				m_vertices[0] = vrt1;
				m_vertices[1] = vrt2;
			}
};


////////////////////////////////////////////////////////////////////////////////////////////////
//	FaceVertices
///	holds the vertices of a Face or an FaceDescriptor
/*
class FaceVertices
{
	friend class Grid;
	public:
		inline VertexBase* vertex(uint index) const	{return m_vertices[index];}
		inline uint num_vertices() const			{return m_vertices.size();}

	//	compatibility with std::vector for some template routines
	///	returns the number of vertices.
		inline size_t size() const	{return m_vertices.size();}
	///	returns the i-th vertex.
		VertexBase* operator[](uint index) const {return m_vertices[index];}

	protected:
		inline void set_num_vertices(int numVrts)	{m_vertices.resize(numVrts);}

	protected:
		typedef std::vector<VertexBase*> 	VertexVec;

	protected:
		VertexVec		m_vertices;
};
*/

//	the following code-section was part of a quick test
//	regarding speed limitations of the originally used std::vector.
//	in a first test the speed increase was quite small.
const int MAX_FACE_VERTICES = 4;
/**	Please note that this class does not have a virtual destructor.*/
class FaceVertices
{
	friend class Grid;
	public:
		inline VertexBase* vertex(uint index) const	{return m_vertices[index];}
		inline size_t num_vertices() const			{return m_numVrts;}

	//	compatibility with std::vector for some template routines
	///	returns the number of vertices.
		inline size_t size() const	{return m_numVrts;}
	///	returns the i-th vertex.
		VertexBase* operator[](size_t index) const {return m_vertices[index];}

	protected:
		inline void set_num_vertices(size_t numVrts)
		{
			assert((numVrts >= 0 && (int)numVrts <= MAX_FACE_VERTICES) && "unsupported number of vertices.");
			m_numVrts = numVrts;
		}

		inline void assign_face_vertices(const FaceVertices& fv)
		{
			m_numVrts = fv.num_vertices();
			for(size_t i = 0; i < m_numVrts; ++i)
				m_vertices[i] = fv.m_vertices[i];
		}

	protected:
		VertexBase* m_vertices[MAX_FACE_VERTICES];
		uint m_numVrts;
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	Face
///	Faces are 2-dimensional objects.
/**
 * Base class for all 2-dimensional objects.
 * Faces connect three or more vertices.
 * A face should always be flat (if viewed in 3 dimensions).
 * You can not create an instance of Face. grids are constructed from derivatives of face.
 * The vertices of a face have always to be specified in counterclockwise order!
 *
 * \ingroup lib_grid_geometric_objects
 */
class Face : public GeometricObject, public FaceVertices
{
	friend class Grid;
	public:
		// lower dimensional Base Object
		typedef EdgeBase lower_dim_base_object;

	public:
		inline static bool type_match(GeometricObject* pObj)	{return dynamic_cast<Face*>(pObj) != NULL;}

		virtual ~Face()	{}

		inline EdgeDescriptor edge(int index) const
			{return EdgeDescriptor(m_vertices[index], m_vertices[(index+1) % size()]);}

		inline void edge(int index, EdgeDescriptor& edOut)
			{edOut.set_vertices(m_vertices[index], m_vertices[(index+1) % size()]);}

		inline uint num_edges() const	{return num_vertices();}
		inline uint num_sides() const	{return num_edges();}

		virtual int shared_pipe_section() const	{return -1;}
		virtual int base_object_type_id() const	{return FACE;}
		virtual int reference_object_id() const	{return -1;}

	/**	A default implementation is featured to allow empty instances of
	 *	this class. This is required to allow the use of this class
	 *	for compile-time method selection by dummy-parameters.
	 *	It is cruical that derived classes overload this method.*/
		virtual EdgeBase* create_edge(int index)	{return NULL;}	///< create the edge with index i and return it.


	/**
	 * The refine method can be used to create new elements by inserting new vertices
	 * on the face.
	 * The user that calls this function is responsible to either register the new
	 * faces with a grid (the grid from which the vertices are), or to take responsibility
	 * for deletion of the acquired memory (delete each element in vNewFacesOut).
	 * - Specify vertices that shall be inserted on edges with newEdgeVertices. Vertices
	 * are inserted on the edge that corresponds to their index. Use NULL to indicate
	 * that no vertex shall be inserted on the associated edge. newEdgeVertices has to point
	 * to an array that holds as many vertices as there are edges in the face.
	 * - If the method has to create a new inner vertex, it will be returned in newFaceVertexOut.
	 * - If you specify pvSubstituteVertices, the created faces will reference the vertices in
	 * pvSubstituteVertices. Note that pvSubstituteVertices has to contain exactly as many
	 * vertices as the refined Face. Vertices with the same index correlate.
	 */
		virtual bool refine(std::vector<Face*>& vNewFacesOut,
							VertexBase** newFaceVertexOut,
							VertexBase** newEdgeVertices,
							VertexBase* newFaceVertex = NULL,
							VertexBase** pSubstituteVertices = NULL)	{return false;}


		/**
		 * The refine_regular method can be used to create new elements by inserting new vertices
		 * on the face.
		 * The user that calls this function is responsible to either register the new
		 * faces with a grid (the grid from which the vertices are), or to take responsibility
		 * for deletion of the acquired memory (delete each element in vNewFacesOut).
		 * - If the user did not specify a newFaceVertex but a face-vertex is required for
		 * regular refinement, a new one is created and a pointer is returned in newFaceVertexOut.
		 * Be sure to either register this vertex at the grid or to delete it.
		 * If no new vertex was created, newFaceVertexOut contains NULL.
		 * - Specify vertices that shall be inserted on edges with vNewEdgeVertices. Vertices
		 * are inserted on the edge that corresponds to their index. For each edge a new
		 * vertex has to be specified.
		 * - If newFaceVertex is not NULL, vNewFaceIndex will be inserted in the center of the face.
		 * - If no newFaceVertex had been specified but a inner vertex is required for regular
		 * refinement, a vertex of the same type as prototypeVertex will be created.
		 * - If you specify pvSubstituteVertices, the created faces will reference the vertices in
		 * pvSubstituteVertices. Note that pvSubstituteVertices has to contain exactly as many
		 * vertices as the refined Face. Vertices with the same index correlate.
		 */
			virtual bool refine_regular(std::vector<Face*>& vNewFacesOut,
								VertexBase** newFaceVertexOut,
								std::vector<VertexBase*>& vNewEdgeVertices,
								VertexBase* newFaceVertex,
								const VertexBase& prototypeVertex,
								VertexBase** pSubstituteVertices = NULL)	{return false;}


	/**
	 * The collapse_edge method creates new geometric objects by collapsing the specified edge.
	 * The user that calls this function is responsible to either register the new
	 * faces with a grid (the grid from which the vertices are), or to take responsibility
	 * for deletion of the acquired memory (delete each element in vNewFacesOut).
	 * - edgeIndex specifies the edge which shall be collapsed. edgeIndex has to lie
	 *   between 0 and this->num_edges().
	 * - Vertices adjacent to the collapsed edge will be replaced by newVertex.
	 * - If you specify pvSubstituteVertices, the created faces will reference the vertices in
	 * pvSubstituteVertices. Note that pvSubstituteVertices has to contain exactly as many
	 * vertices as the face in which you collapse the edge. Vertices with the same index correlate.
	 */
		virtual bool collapse_edge(std::vector<Face*>& vNewFacesOut,
								int edgeIndex, VertexBase* newVertex,
								VertexBase** pSubstituteVertices = NULL)	{return false;}

	/**
	 * The collapse_edgea method creates new geometric objects by collapsing the specified edges
	 * simultaneously. This method makes sense only for faces with more than 4 edges.
	 * The user that calls this function is responsible to either register the new
	 * faces with a grid (the grid from which the vertices are), or to take responsibility
	 * for deletion of the acquired memory (delete each element in vNewFacesOut).
	 * - for each entry in vNewEdgeVertices which is not NULL, the edge with the same index will
	 *    be collapsed and replaced by the specified vertex.
	 *    The size of vNewEdgeVertices has thus to be between 0 and this->num_edges().
	 * - If you specify pvSubstituteVertices, the created faces will reference the vertices in
	 * pvSubstituteVertices. Note that pvSubstituteVertices has to contain exactly as many
	 * vertices as the face in which you collapse the edge. Vertices with the same index correlate.
	 */
		virtual bool collapse_edges(std::vector<Face*>& vNewFacesOut,
								std::vector<VertexBase*>& vNewEdgeVertices,
								VertexBase** pSubstituteVertices = NULL)	{return false;}

// BEGIN Depreciated
	/**	creates the faces that result from the splitting of the edge with index 'splitEdgeIndex'.
	 *  The user that calls this function is responsible to either register those faces with the
	 *  same grid in which 'this' was registered, or to take over responsibility for deletion of the
	 *  acquired memory.
	 *  With pvSubstituteVertices you can optionally pass a vector of vertices, which will be used
	 *  instead of the original vertices of this face. If specified, pvSubstituteVertices has to
	 *  contain exactly as many vertices as the face itself.
	 */
		virtual void create_faces_by_edge_split(int splitEdgeIndex,
							VertexBase* newVertex,
							std::vector<Face*>& vNewFacesOut,
							VertexBase** pSubstituteVertices = NULL)	{};
// END Depreciated

	/**	creates the faces that result from the collapsing of the edge with index 'splitEdgeIndex'.*/
		//virtual void create_faces_by_edge_collapse(int collapseEdgeIndex,
		//						VertexBase* newVertex,
		//						std::vector<Face*>& vNewFacesOut) = 0;

	protected:
		inline void set_vertex(uint index, VertexBase* pVrt)	{m_vertices[index] = pVrt;}
};

template <>
class geometry_traits<Face>
{
	public:
		typedef GenericGeometricObjectIterator<Face*>		iterator;
		typedef ConstGenericGeometricObjectIterator<Face*>	const_iterator;

		typedef Face	geometric_base_object;
		//typedef void Descriptor;	///< Faces can't be created directly

		enum
		{
			SHARED_PIPE_SECTION = -1,
			BASE_OBJECT_TYPE_ID = FACE
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_INVALID;
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	FaceDescriptor
///	Can be queried for the edges and vertices of a face.
class FaceDescriptor : public FaceVertices
{
	public:
		FaceDescriptor();
		FaceDescriptor(uint numVertices);
		FaceDescriptor(const FaceDescriptor& fd);

		FaceDescriptor& operator = (const FaceDescriptor& fd);

		inline void set_num_vertices(uint numVertices)	{FaceVertices::set_num_vertices(numVertices);}
		inline void set_vertex(uint index, VertexBase* vrt)
			{m_vertices[index] = vrt;}

/*
		inline EdgeDescriptor edge(int index) const
			{return EdgeDescriptor(m_vertices[index], m_vertices[(index+1) % m_vertices.size()]);}

		inline void edge(int index, EdgeDescriptor& edOut) const
			{edOut.set_vertices(m_vertices[index], m_vertices[(index+1) % m_vertices.size()]);}

		inline uint num_edges() const	{return m_vertices.size();}
*/
};


////////////////////////////////////////////////////////////////////////////////////////////////
//	VolumeVertices
///	holds the vertices of a Volume or a VolumeDescriptor
/**	Please note that this class does not have a virtual destructor.*/
class VolumeVertices
{
	friend class Grid;
	public:
		inline VertexBase* vertex(uint index) const	{return m_vertices[index];}
		inline uint num_vertices() const			{return m_vertices.size();}

	//	compatibility with std::vector for some template routines
	///	returns the number of vertices.
		inline size_t size() const	{return m_vertices.size();}
	///	returns the i-th vertex.
		VertexBase* operator[](uint index) const {return m_vertices[index];}

	protected:
		inline void assign_volume_vertices(const VolumeVertices& vv)
		{
			m_vertices.resize(vv.num_vertices());
			for(size_t i = 0; i < num_vertices(); ++i)
				m_vertices[i] = vv.m_vertices[i];
		}

	protected:
		typedef std::vector<VertexBase*> 	VertexVec;

	protected:
		VertexVec		m_vertices;
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	Volume
///	Volumes are 3-dimensional objects.
/**
 * Base class for all 3-dimensional objects.
 * Volumes connect four or more vertices.
 *
 * default implementations of all methods are featured to allow
 * empty instances of this class.
 * This is required to allow the use of this class
 * for compile-time method selection by dummy-parameters.
 * It is cruical that derived classes overload thoes methods.
 *
 * \ingroup lib_grid_geometric_objects
 */
class Volume : public GeometricObject, public VolumeVertices
{
	friend class Grid;
	public:
		// lower dimensional Base Object
		typedef Face lower_dim_base_object;

	public:
		inline static bool type_match(GeometricObject* pObj)	{return dynamic_cast<Volume*>(pObj) != NULL;}

		virtual ~Volume()	{}

		virtual EdgeDescriptor edge(int index) const				{return EdgeDescriptor(NULL, NULL);}
		virtual void edge(int index, EdgeDescriptor& edOut) const	{edOut = EdgeDescriptor(NULL, NULL);}
		virtual uint num_edges() const								{return 0;}

		virtual FaceDescriptor face(int index) const				{return FaceDescriptor(0);}
		virtual void face(int index, FaceDescriptor& fdOut) const	{fdOut = FaceDescriptor(0);}
		virtual uint num_faces() const								{return 0;}
		inline uint num_sides() const								{return num_faces();}

	/**
	 * The refine method can be used to create new elements by inserting new vertices
	 * on the volume.
	 *
	 * New volumes will be returned in vNewFacesOut.
	 * If a new vertex has to be created from the prototypeVertex (this happens in
	 * more complicated situations) the pointer to the new vertex is returned in
	 * ppNewVertexOut. ppNewVertexOut contains NULL if no new vertex has been created.
	 *
	 * The user that calls this function is responsible to either register the new
	 * volumes and the new vertex with a grid (the grid from which the referenced
	 * vertices are), or to take responsibility for deletion of the acquired memory
	 * (delete each element in vNewVolumesOut and ppNewVertexOut - if it is not NULL).
	 *
	 * - Specify vertices that shall be inserted on edges with vNewEdgeVertices. Vertices
	 *   are inserted on the edge that corresponds to its index. Use NULL to indicate
	 *   that no vertex shall be inserted on the associated edge.
	 * - Specify vertices that shall be inserted on faces with vNewFaceVertices. Vertices
	 *   are inserted on the face that corresponds to its index. Use NULL to indicate
	 *   that no vertex shall be inserted on the associated face.
	 * - If newVolumeVertex is not NULL, newVolumeVertex will be inserted in the center of
	 *   the volume.
	 * - via the prototypeVertex you can specify the vertex-type of the vertex that is
	 *   auto-inserted if the refine operation is too complex. In most cases this will be
	 *   an instance of a standard vertex (a concrete type is required. VertexBase will not do).
	 *   The prototypeVertex has not to be registered at any grid - it may be a temporary
	 *   instance.
	 * - If you specify pvSubstituteVertices, the created volumes will reference the vertices
	 *   in pvSubstituteVertices. Note that pvSubstituteVertices has to contain exactly as
	 *   many vertices as the refined Face. Vertices with the same index correlate.
	 */
		virtual bool refine(std::vector<Volume*>& vNewVolumesOut,
							VertexBase** ppNewVertexOut,
							VertexBase** newEdgeVertices,
							VertexBase** newFaceVertices,
							VertexBase* newVolumeVertex,
							const VertexBase& prototypeVertex,
							VertexBase** pSubstituteVertices = NULL)	{return false;}
	/**
	 * See the refine-method for a detailed explanation of each parameter.
	 * Please note that for each edge and for each face a new vertex has to be specified.
	 * This means that vNewEdgeVertices and vNewFaceVertices may not contain NULL-pointers.
	 */
		virtual bool refine_regular(std::vector<Volume*>& vNewVolumesOut,
							VertexBase** ppNewVertexOut,
							std::vector<VertexBase*>& vNewEdgeVertices,
							std::vector<VertexBase*>& vNewFaceVertices,
							VertexBase* newVolumeVertex,
							const VertexBase& prototypeVertex,
							std::vector<VertexBase*>* pvSubstituteVertices = NULL)	{return false;}

	/**
	 * The collapse_edge method creates new geometric objects by collapsing the specified edge.
	 * The user that calls this function is responsible to either register the new
	 * faces with a grid (the grid from which the vertices are), or to take responsibility
	 * for deletion of the acquired memory (delete each element in vNewFacesOut).
	 * - edgeIndex specifies the edge which shall be collapsed. edgeIndex has to lie
	 *   between 0 and this->num_edges().
	 * - Vertices adjacent to the collapsed edge will be replaced by newVertex.
	 * - If you specify pvSubstituteVertices, the created faces will reference the vertices in
	 * pvSubstituteVertices. Note that pvSubstituteVertices has to contain exactly as many
	 * vertices as the face in which you collapse the edge. Vertices with the same index correlate.
	 */
		virtual bool collapse_edge(std::vector<Volume*>& vNewVolumesOut,
								int edgeIndex, VertexBase* newVertex,
								std::vector<VertexBase*>* pvSubstituteVertices = NULL)	{return false;}

		virtual int shared_pipe_section() const	{return -1;}
		virtual int base_object_type_id() const	{return VOLUME;}
		virtual int reference_object_id() const	{return -1;}

	protected:
		virtual EdgeBase* create_edge(int index)	{return NULL;}	///< create the edge with index i and return it.
		virtual Face* create_face(int index)		{return NULL;}	///< create the face with index i and return it.

	/**	creates the volumes that result from the splitting of the edge with index 'splitEdgeIndex'.*/
		//virtual void create_volumes_by_edge_split(int splitEdgeIndex,
		//						VertexBase* newVertex,
		//						std::vector<Volume*>& vNewFacesOut) = 0;

	/**	creates the volumes that result from the collapsing of the edge with index 'splitEdgeIndex'.*/
		//virtual void create_Volumes_by_edge_collapse(int collapseEdgeIndex,
		//						VertexBase* newVertex,
		//						std::vector<Volume*>& vNewFacesOut) = 0;

		inline void set_vertex(uint index, VertexBase* pVrt)	{m_vertices[index] = pVrt;}
};

template <>
class geometry_traits<Volume>
{
	public:
		typedef GenericGeometricObjectIterator<Volume*>			iterator;
		typedef ConstGenericGeometricObjectIterator<Volume*>	const_iterator;

		typedef Volume		geometric_base_object;

		enum
		{
			SHARED_PIPE_SECTION = -1,
			BASE_OBJECT_TYPE_ID = VOLUME
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_INVALID;
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	VolumeDescriptor
///	Can be queried for the edges, faces and vertices of a volume.
class VolumeDescriptor : public VolumeVertices
{
	public:
		VolumeDescriptor();
		VolumeDescriptor(uint numVertices, uint numEdges, uint numFaces);
		VolumeDescriptor(const VolumeDescriptor& vd);

		VolumeDescriptor& operator = (const VolumeDescriptor& vd);

		inline void set_num_vertices(uint numVertices)	{m_vertices.resize(numVertices);}
		inline void set_vertex(uint index, VertexBase* vrt)
			{m_vertices[index] = vrt;}

/*
		inline EdgeDescriptor edge(uint index) const				{return m_edges[index];}
		inline void edge(uint index, EdgeDescriptor& edOut) const	{edOut = m_edges[index];}
		inline uint num_edges() const								{return m_edges.size();}
		inline void set_num_edges(uint numEdges)				{m_edges.resize(numEdges);}
		inline void set_edge_descriptor(uint index, const EdgeDescriptor& ed)
			{m_edges[index] = ed;}

		inline FaceDescriptor face(uint index) const				{return m_faces[index];}
		inline void face(uint index, FaceDescriptor& fdOut)	const	{fdOut = m_faces[index];}
		inline uint num_faces() const								{return m_faces.size();}
		inline void set_num_faces(uint numFaces)				{m_faces.resize(numFaces);}
		inline void set_face_descriptor(uint index, const FaceDescriptor& fd)
			{m_faces[index] = fd;}

	protected:
		typedef std::vector<EdgeDescriptor>	EdgeDescriptorVec;
		typedef std::vector<FaceDescriptor>	FaceDescriptorVec;

	protected:
		EdgeDescriptorVec	m_edges;
		FaceDescriptorVec	m_faces;
*/
};


typedef geometry_traits<VertexBase>::iterator	VertexBaseIterator;
typedef geometry_traits<EdgeBase>::iterator		EdgeBaseIterator;
typedef geometry_traits<Face>::iterator			FaceIterator;
typedef geometry_traits<Volume>::iterator		VolumeIterator;

typedef geometry_traits<VertexBase>::const_iterator	ConstVertexBaseIterator;
typedef geometry_traits<EdgeBase>::const_iterator	ConstEdgeBaseIterator;
typedef geometry_traits<Face>::const_iterator		ConstFaceIterator;
typedef geometry_traits<Volume>::const_iterator		ConstVolumeIterator;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	specialization of attachment_traits for libGrid::GeometricObject

template<>
class attachment_traits<GeometricObject*, Grid>
{
	public:
		typedef GeometricObject*&		ElemRef;
		typedef GeometricObject*		ElemPtr;
		typedef const GeometricObject*	ConstElemPtr;
		typedef Grid*					ElemHandlerPtr;
		typedef const Grid*				ConstElemHandlerPtr;

		static inline void invalidate_entry(ElemHandlerPtr pHandler, ElemRef elem)			{elem = NULL;}
		static inline bool entry_is_invalid(ElemHandlerPtr pHandler, ElemRef elem)			{return elem != NULL;}
		static inline uint get_data_index(ConstElemHandlerPtr pHandler, ConstElemPtr elem)	{return elem->m_gridDataIndex;}
		static inline void set_data_index(ElemHandlerPtr pHandler, ElemPtr elem, uint index){elem->m_gridDataIndex = index;}
};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	hash-funtions for vertices
///	returns the hash-value of the vertex.
template <>
unsigned long hash_key<PVertexBase>(const PVertexBase& key);

////////////////////////////////////////////////////////////////////////
//	hash-funtions for edges
///	the hash-key is a function of vertex-hash-values.
/**
 * The hash value depends on the associated vertices.
 * If an EdgeBase (or EdgeDescriptor) has the same vertices
 * as another EdgeBase (or EdgeDescriptor), the hash-keys
 * are the same.
 */
template <>
unsigned long hash_key<PEdgeVertices>(const PEdgeVertices& key);

///	the hash-key is a function of vertex-hash-values.
/** \sa hash_key<PEdgeVertices>*/
template <>
unsigned long hash_key<PEdgeBase>(const PEdgeBase& key);

///	the hash-key is a function of vertex-hash-values.
/** \sa hash_key<PEdgeVertices>*/
template <>
unsigned long hash_key<PEdgeDescriptor>(const PEdgeDescriptor& key);

////////////////////////////////////////////////////////////////////////
//	hash-funtions for faces
///	the hash-key is a function of vertex-hash-values.
/**
 * The hash value depends on the associated vertices.
 * If an Face (or FaceDescriptor) has the same vertices
 * as another Face (or FaceDescriptor), the hash-keys
 * are the same.
 */
template <>
unsigned long hash_key<PFaceVertices>(const PFaceVertices& key);

///	the hash-key is a function of vertex-hash-values.
/**\sa hash_key<PFaceVertices>*/
template <>
unsigned long hash_key<PFace>(const PFace& key);

///	the hash-key is a function of vertex-hash-values.
/**\sa hash_key<PFaceVertices>*/
template <>
unsigned long hash_key<PFaceDescriptor>(const PFaceDescriptor& key);

////////////////////////////////////////////////////////////////////////
//	hash-funtions for volumes
///	the hash-key is a function of vertex-hash-values.
/**
 * The hash value depends on the associated vertices.
 * If an Volume (or VolumeDescriptor) has the same vertices
 * as another Volume (or VolumeDescriptor), the hash-keys
 * are the same.
 */
template <>
unsigned long hash_key<PVolumeVertices>(const PVolumeVertices& key);

///	the hash-key is a function of vertex-hash-values.
/**\sa hash_key<PVolumeVertices>*/
template <>
unsigned long hash_key<PVolume>(const PVolume& key);

///	the hash-key is a function of vertex-hash-values.
/**\sa hash_key<PVolumeVertices>*/
template <>
unsigned long hash_key<PVolumeDescriptor>(const PVolumeDescriptor& key);

}//	end of namespace

#endif
