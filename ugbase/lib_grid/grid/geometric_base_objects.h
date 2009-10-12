//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d09

#ifndef __H__LIB_GRID__GEOMETRIC_BASE_OBJECTS__
#define __H__LIB_GRID__GEOMETRIC_BASE_OBJECTS__

#include <list>
#include "common/types.h"
#include "common/util/attachment_pipe.h"

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
//	Predeclaration of geometric objects.
class GeometricObject;
class VertexBase;
class EdgeBase;
class Face;
class Volume;

////////////////////////////////////////////////////////////////////////
//	GeometricObjectContainer
///	Declaration of the container that will hold all geometric objects.
typedef std::list<GeometricObject*> 		GeometricObjectContainer;

////////////////////////////////////////////////////////////////////////
//	GeometricObjectIterator
///	This Iterator will be used as base-class for iterators of specialized geometric objects.
typedef GeometricObjectContainer::iterator	GeometricObjectIterator;

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

	public:
		typedef TValue	value_type;

	public:
		GenericGeometricObjectIterator()	{}

		GenericGeometricObjectIterator(const GenericGeometricObjectIterator& iter) :
			TBaseIterator(iter)	{}

		GenericGeometricObjectIterator(const typename std::list<TValue>::iterator& iter) :
			TBaseIterator(*((TBaseIterator*)&iter)){}

		inline TValue& operator* ()	{return (TValue&)(GeometricObjectContainer::iterator::operator*());}

	protected:
		/**
		 * This method should only be invoked by Grid.
		 * If you get an error in your code pointing to this method, then you are
		 * most presumably assigning an iterator of the wrong type.
		 */
		/*
		GenericGeometricObjectIterator(const GeometricObjectContainer::iterator& iter) :
			GeometricObjectIterator(iter){}
		*/
};

////////////////////////////////////////////////////////////////////////
//	iterator_cast
///	This method is unsave and should not be used until you know exactly what you are doing!
template <class TIterDest, class TIterSrc>
inline TIterDest
iterator_cast(const TIterSrc& iter)
{
	return *reinterpret_cast<const TIterDest*>(&iter);
}

/**
 * \defgroup GeometricObjects Geometric Objects
 * \brief Geometric objects are the building blocks of a grid.
 */
////////////////////////////////////////////////////////////////////////////////////////////////
//	GeometricObject
///	The base class for all geometric objects, such as vertices, edges, faces, volumes, ...
/**
 * In order to be used by libGrid, all derivatives of GeometricObject
 * have to specialize geometry_traits<GeomObjectType>.
 *
 * \ingroup GeometricObjects
 */
class GeometricObject
{
	friend class Grid;
	friend class attachment_traits<GeometricObject*, Grid>;
	public:
		virtual ~GeometricObject()	{}

	///	create an instance of the derived type
		virtual GeometricObject* create_empty_instance() const = 0;

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
		typedef GeometricObjectIterator		iterator;
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
 * \ingroup GeometricObjects
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
};

template <>
class geometry_traits<VertexBase>
{
	public:
		typedef GenericGeometricObjectIterator<VertexBase*>	iterator;
		typedef VertexBase	geometric_base_object;

		enum
		{
			SHARED_PIPE_SECTION = -1,
			BASE_OBJECT_TYPE_ID = VERTEX
		};
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	EdgeBase
///	Base-class for edges
/**
 * EdgeBase is the base class of all 1-dimensional geometric objects.
 * Edges connect two vertices.
 *
 * \ingroup GeometricObjects
 */
class EdgeBase : public GeometricObject
{
	friend class Grid;
	public:
		inline static bool type_match(GeometricObject* pObj)	{return dynamic_cast<EdgeBase*>(pObj) != NULL;}

		virtual ~EdgeBase()	{}

		inline VertexBase* vertex(uint index)	{return m_vertices[index];}
		inline uint num_vertices()				{return 2;}	// this method is supplied to allow the use of EdgeBase in template-methods that require a num_vertices() method.

		virtual int shared_pipe_section() const	{return -1;}
		virtual int base_object_type_id() const	{return EDGE;}
		virtual int reference_object_id() const	{return -1;}

	/**
	 * create 2 new edges, connecting the original edges end-points with vrtNew.
	 * Newly created edges have to be registered at a grid manually by the caller.
	 * If the caller does not register the edges in vGeomOut at a grid, he is
	 * responsible to free the associated memory (delete each element in vNewEdgesOut).
	 * Please note that refining an edge using this method does not automatically
	 * refine associated elements.
	 * Be sure to store the new edges in the right order. vNewEdgesOut should contain
	 * the edge connecting vertex(0) and newVertex first.
	 */
		virtual bool refine(std::vector<EdgeBase*>& vNewEdgesOut,
											VertexBase* newVertex)	{return false;}

	protected:
		inline void set_vertex(uint index, VertexBase* pVrt)	{m_vertices[index] = pVrt;}
		VertexBase*	m_vertices[2];
};

template <>
class geometry_traits<EdgeBase>
{
	public:
		typedef GenericGeometricObjectIterator<EdgeBase*>	iterator;
		typedef EdgeBase	geometric_base_object;

		enum
		{
			SHARED_PIPE_SECTION = -1,
			BASE_OBJECT_TYPE_ID = EDGE
		};
};


////////////////////////////////////////////////////////////////////////////////////////////////
//	EdgeDescriptor
///	Can be used to store information about an edge and to construct an edge.
class EdgeDescriptor
{
	public:
		EdgeDescriptor()	{}
		EdgeDescriptor(const EdgeDescriptor& ed)
			{
				m_vertex[0] = ed.vertex(0);
				m_vertex[1] = ed.vertex(1);
			}
		EdgeDescriptor(VertexBase* vrt1, VertexBase* vrt2)
			{
				m_vertex[0] = vrt1;
				m_vertex[1] = vrt2;
			}

		inline VertexBase* vertex(uint index) const	{return m_vertex[index];}
		inline void set_vertex(uint index, VertexBase* vrt)	{m_vertex[index] = vrt;}
		inline void set_vertices(VertexBase* vrt1, VertexBase* vrt2)
			{
				m_vertex[0] = vrt1;
				m_vertex[1] = vrt2;
			}

	protected:
		VertexBase* m_vertex[2];
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
 * \ingroup GeometricObjects
 */
class Face : public GeometricObject
{
	friend class Grid;
	public:
		inline static bool type_match(GeometricObject* pObj)	{return dynamic_cast<Face*>(pObj) != NULL;}

		virtual ~Face()	{}

		inline VertexBase* vertex(uint index)	{return m_vertices[index];}
		inline uint num_vertices()			{return m_vertices.size();}

		inline EdgeDescriptor edge(int index) const
			{return EdgeDescriptor(m_vertices[index], m_vertices[(index+1) % m_vertices.size()]);}

		inline void edge(int index, EdgeDescriptor& edOut)
			{edOut.set_vertices(m_vertices[index], m_vertices[(index+1) % m_vertices.size()]);}

		inline uint num_edges()	{return m_vertices.size();}

		virtual int shared_pipe_section() const	{return -1;}
		virtual int base_object_type_id() const	{return FACE;}
		virtual int reference_object_id() const	{return -1;}

		virtual EdgeBase* create_edge(int index) = 0;	///< create the edge with index i and return it.


	/**
	 * The refine method can be used to create new elements by inserting new vertices
	 * on the face.
	 * The user that calls this function is responsible to either register the new
	 * faces with a grid (the grid from which the vertices are), or to take responsibility
	 * for deletion of the acquired memory (delete each element in vNewFacesOut).
	 * - Specify vertices that shall be inserted on edges with vNewEdgeVertices. Vertices
	 * are inserted on the edge that corresponds to their index. Use NULL to indicate
	 * that no vertex shall be inserted on the associated edge.
	 * - If newFaceVertex is not NULL, vNewFaceIndex will be inserted in the center of the face.
	 * - If you specify pvSubstituteVertices, the created faces will reference the vertices in
	 * pvSubstituteVertices. Note that pvSubstituteVertices has to contain exactly as many
	 * vertices as the refined Face. Vertices with the same index correlate.
	 */
		virtual bool refine(std::vector<Face*>& vNewFacesOut,
							std::vector<VertexBase*>& vNewEdgeVertices,
							VertexBase* newFaceVertex,
							std::vector<VertexBase*>* pvSubstituteVertices = NULL)	{return false;}


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
		virtual bool collapse_edge(std::vector<Face*>& vNewFacesOut,
								int edgeIndex, VertexBase* newVertex,
								std::vector<VertexBase*>* pvSubstituteVertices = NULL)	{return false;}

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
								std::vector<VertexBase*>* pvSubstituteVertices = NULL)	{return false;}

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
							std::vector<VertexBase*>* pvSubstituteVertices = NULL) = 0;
// END Depreciated

	/**	creates the faces that result from the collapsing of the edge with index 'splitEdgeIndex'.*/
		//virtual void create_faces_by_edge_collapse(int collapseEdgeIndex,
		//						VertexBase* newVertex,
		//						std::vector<Face*>& vNewFacesOut) = 0;

	protected:
		inline void set_vertex(uint index, VertexBase* pVrt)	{m_vertices[index] = pVrt;}

	protected:
		typedef std::vector<VertexBase*> 	VertexVec;

	protected:
		VertexVec		m_vertices;
};

template <>
class geometry_traits<Face>
{
	public:
		typedef GenericGeometricObjectIterator<Face*>	iterator;
		typedef Face	geometric_base_object;
		//typedef void Descriptor;	///< Faces can't be created directly

		enum
		{
			SHARED_PIPE_SECTION = -1,
			BASE_OBJECT_TYPE_ID = FACE
		};
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	FaceDescriptor
///	Can be queried for the edges and vertices of a face.
class FaceDescriptor
{
	public:
		FaceDescriptor()	{}

		FaceDescriptor(uint numVertices)
			{
				set_num_vertices(numVertices);
			}

		virtual ~FaceDescriptor()	{}

		inline VertexBase* vertex(uint index) const		{return m_vertices[index];}
		inline uint num_vertices() const				{return m_vertices.size();}
		inline void set_num_vertices(uint numVertices)	{m_vertices.resize(numVertices);}
		inline void set_vertex(uint index, VertexBase* vrt)
			{m_vertices[index] = vrt;}


		inline EdgeDescriptor edge(int index) const
			{return EdgeDescriptor(m_vertices[index], m_vertices[(index+1) % m_vertices.size()]);}

		inline void edge(int index, EdgeDescriptor& edOut) const
			{edOut.set_vertices(m_vertices[index], m_vertices[(index+1) % m_vertices.size()]);}

		inline uint num_edges() const	{return m_vertices.size();}

	protected:
		typedef std::vector<VertexBase*> 	VertexVec;

	protected:
		VertexVec			m_vertices;
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	Volume
///	Volumes are 3-dimensional objects.
/**
 * Base class for all 3-dimensional objects.
 * Volumes connect four or more vertices.
 *
 * \ingroup GeometricObjects
 */
class Volume : public GeometricObject
{
	friend class Grid;
	public:
		inline static bool type_match(GeometricObject* pObj)	{return dynamic_cast<Volume*>(pObj) != NULL;}

		virtual ~Volume()	{}

		inline VertexBase* vertex(uint index)	{return m_vertices[index];}
		inline uint num_vertices()				{return m_vertices.size();}

		virtual EdgeDescriptor edge(int index) const = 0;
		virtual void edge(int index, EdgeDescriptor& edOut) const = 0;
		virtual uint num_edges() const = 0;

		virtual FaceDescriptor face(int index) const  = 0;
		virtual void face(int index, FaceDescriptor& fdOut) const = 0;
		virtual uint num_faces() const = 0;

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
							std::vector<VertexBase*>& vNewEdgeVertices,
							std::vector<VertexBase*>& vNewFaceVertices,
							VertexBase* newVolumeVertex,
							const VertexBase& prototypeVertex,
							std::vector<VertexBase*>* pvSubstituteVertices = NULL)	{return false;}

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
		virtual EdgeBase* create_edge(int index) = 0;	///< create the edge with index i and return it.
		virtual Face* create_face(int index) = 0;		///< create the face with index i and return it.

	/**	creates the volumes that result from the splitting of the edge with index 'splitEdgeIndex'.*/
		//virtual void create_volumes_by_edge_split(int splitEdgeIndex,
		//						VertexBase* newVertex,
		//						std::vector<Volume*>& vNewFacesOut) = 0;

	/**	creates the volumes that result from the collapsing of the edge with index 'splitEdgeIndex'.*/
		//virtual void create_Volumes_by_edge_collapse(int collapseEdgeIndex,
		//						VertexBase* newVertex,
		//						std::vector<Volume*>& vNewFacesOut) = 0;

		inline void set_vertex(uint index, VertexBase* pVrt)	{m_vertices[index] = pVrt;}

	protected:
		typedef std::vector<VertexBase*> 	VertexVec;

	protected:
		VertexVec		m_vertices;
};

template <>
class geometry_traits<Volume>
{
	public:
		typedef GenericGeometricObjectIterator<Volume*>		iterator;
		typedef Volume		geometric_base_object;

		enum
		{
			SHARED_PIPE_SECTION = -1,
			BASE_OBJECT_TYPE_ID = VOLUME
		};
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	VolumeDescriptor
///	Can be queried for the edges, faces and vertices of a volume.
class VolumeDescriptor
{
	public:
		VolumeDescriptor()	{}

		VolumeDescriptor(uint numVertices, uint numEdges, uint numFaces)
			{
				set_num_vertices(numVertices);
				set_num_edges(numEdges);
				set_num_faces(numFaces);
			}

		virtual ~VolumeDescriptor()	{}

		inline VertexBase* vertex(uint index) const		{return m_vertices[index];}
		inline uint num_vertices() const				{return m_vertices.size();}
		inline void set_num_vertices(uint numVertices)	{m_vertices.resize(numVertices);}
		inline void set_vertex(uint index, VertexBase* vrt)
			{m_vertices[index] = vrt;}


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
		typedef std::vector<VertexBase*> 	VertexVec;
		typedef std::vector<EdgeDescriptor>	EdgeDescriptorVec;
		typedef std::vector<FaceDescriptor>	FaceDescriptorVec;

	protected:
		VertexVec			m_vertices;
		EdgeDescriptorVec	m_edges;
		FaceDescriptorVec	m_faces;
};


typedef geometry_traits<VertexBase>::iterator	VertexBaseIterator;
typedef geometry_traits<EdgeBase>::iterator		EdgeBaseIterator;
typedef geometry_traits<Face>::iterator			FaceIterator;
typedef geometry_traits<Volume>::iterator		VolumeIterator;


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	specialization of attachment_traits for libGrid::GeometricObject

template<>
class attachment_traits<GeometricObject*, Grid>
{
	public:
		static inline uint get_data_index(const GeometricObject* elem)					{return elem->m_gridDataIndex;}
		static inline void set_data_index(GeometricObject* elem, uint index)			{elem->m_gridDataIndex = index;};
};

}//	end of namespace

#endif
