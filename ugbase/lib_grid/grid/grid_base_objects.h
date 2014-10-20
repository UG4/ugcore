//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d09

#ifndef __H__LIB_GRID__GEOMETRIC_BASE_OBJECTS__
#define __H__LIB_GRID__GEOMETRIC_BASE_OBJECTS__

#include <list>
#include <cassert>
#include <iostream>
#include <utility>
#include "common/types.h"
#include "common/common.h"
#include "common/assert.h"
#include "lib_grid/attachments/attachment_pipe.h"
#include "lib_grid/attachments/attached_list.h"
#include "common/util/hash_function.h"
#include "common/allocators/small_object_allocator.h"
#include "common/math/ugmath_types.h"
#include "common/util/pointer_const_array.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	CONSTANTS

////////////////////////////////////////////////////////////////////////
///	enumeration of the GridBaseObjects that make up a grid.
enum GridBaseObjectId
{
	VERTEX = 0,
	EDGE,
	FACE,
	VOLUME,
	NUM_GEOMETRIC_BASE_OBJECTS		//always last!!!
};

////////////////////////////////////////////////////////////////////////
//	Reference-Object IDs
///	these ids are used to identify the shape of a geometric object.
enum ReferenceObjectID
{
	ROID_UNKNOWN = -1,
	ROID_VERTEX,
	ROID_EDGE,
	ROID_TRIANGLE,
	ROID_QUADRILATERAL,
	ROID_TETRAHEDRON,
	ROID_HEXAHEDRON,
	ROID_PRISM,
	ROID_PYRAMID,
	ROID_OCTAHEDRON,
	NUM_REFERENCE_OBJECTS
};

inline
ReferenceObjectID operator++(ReferenceObjectID& roid, int)
{
	int tmp = roid;
	return roid = static_cast<ReferenceObjectID>(++tmp);
};

////////////////////////////////////////////////////////////////////////
inline
std::ostream& operator<< (std::ostream& outStream, ReferenceObjectID type)
{
	switch(type)
	{
		case ROID_UNKNOWN: outStream << "(invalid)"; break;
		case ROID_VERTEX: outStream << "Vertex"; break;
		case ROID_EDGE: outStream << "Edge"; break;
		case ROID_TRIANGLE: outStream << "Triangle"; break;
		case ROID_QUADRILATERAL: outStream << "Quadrilateral"; break;
		case ROID_TETRAHEDRON: outStream << "Tetrahedron"; break;
		case ROID_OCTAHEDRON: outStream << "Octahedron"; break;
		case ROID_HEXAHEDRON: outStream << "Hexahedron"; break;
		case ROID_PRISM: outStream << "Prism"; break;
		case ROID_PYRAMID: outStream << "Pyramid"; break;
		default: throw(UGError("Unknown ReferenceObjectID in operator<<"));
	}
	return outStream;
};


////////////////////////////////////////////////////////////////////////
//	PREDECLARATIONS
class Grid;
template<class TElem>	class ElementStorage;

class GridObject;	//	geometric base object
class Vertex;		//	base for all 0-dimensional grid objects.
class Edge;			//	base for all 1-dimensional grid objects.
class Face;				//	base for all 2-dimensional grid objects.
class Volume;			//	base for all 3-dimensional grid objects.

class EdgeDescriptor;	//	describes an edge.
class FaceDescriptor;	//	describes a face.
class VolumeDescriptor;	//	describes a volume.

class EdgeVertices;		//	manages the vertices of an edge. Base for Edge and EdgeDescriptor.
class FaceVertices;		//	manages the vertices of a face. Base for Face and FaceDescriptor.
class VolumeVertices;	//	manages the vertices of a volume. Base for Volume and VolumeDescriptor.

//	pointer-types. Primarily required for template-specializations.
typedef Vertex*	PVertex;
typedef Edge*		PEdge;
typedef Face*			PFace;
typedef Volume*		PVolume;

typedef EdgeDescriptor*	PEdgeDescriptor;
typedef FaceDescriptor*	PFaceDescriptor;
typedef VolumeDescriptor*	PVolumeDescriptor;

typedef EdgeVertices*		PEdgeVertices;
typedef FaceVertices*		PFaceVertices;
typedef VolumeVertices*		PVolumeVertices;

template<> class attachment_traits<Vertex*, ElementStorage<Vertex> >;
template<> class attachment_traits<Edge*, ElementStorage<Edge> >;
template<> class attachment_traits<Face*, ElementStorage<Face> >;
template<> class attachment_traits<Volume*, ElementStorage<Volume> >;

/**
 * \brief Geometric objects are the building blocks of a grid.
 *
 * \defgroup lib_grid_grid_objects geometric objects
 * \ingroup lib_grid
 */
////////////////////////////////////////////////////////////////////////
//	GridObject
///	The base class for all geometric objects, such as vertices, edges, faces, volumes, ...
/**
 * In order to be used by libGrid, all derivatives of GridObject
 * have to specialize geometry_traits<GeomObjectType>.
 *
 * \ingroup lib_grid_grid_objects
 */
class UG_API GridObject/* : public SmallObject<>*/
{
	friend class Grid;
	friend class attachment_traits<Vertex*, ElementStorage<Vertex> >;
	friend class attachment_traits<Edge*, ElementStorage<Edge> >;
	friend class attachment_traits<Face*, ElementStorage<Face> >;
	friend class attachment_traits<Volume*, ElementStorage<Volume> >;

	public:
		virtual ~GridObject()	{}

	///	create an instance of the derived type
	/**	Make sure to overload this method in derivates of this class!*/
		virtual GridObject* create_empty_instance() const {return NULL;}

		virtual int container_section() const = 0;
		virtual int base_object_id() const = 0;
	/**
	 * A reference object represents a class of geometric objects.
	 * Tetrahedrons, Triangles etc are such classes.
	 * Reference ids should be defined in the file in which concrete geometric objects are defined.
	 */
		virtual ReferenceObjectID reference_object_id() const = 0;///	returns the id of the reference-object.

	///	returns true if the object constrains other objects.
	/**	This is normally only the case for special constraining objects.
	 * The default implementation returns false.*/
		virtual bool is_constraining() const					{return false;}

	///	returns true if the object is constrained by other objects.
	/**	This is normally only the case for special constrained objects.
	 * The default implementation returns false.*/
		virtual bool is_constrained() const						{return false;}

	protected:
	///	ATTENTION: Use this method with extreme care!
	/**	This method is for internal use only and should almost never be called
	 * by a user of lib_grid. The method sets the attachment data index and is
	 * mainly used by attachment-traits classes.*/
		inline void set_grid_data_index(uint index)		{m_gridDataIndex = index;}

	///	Returns the grid attachment data index of a geometric object.
		inline uint grid_data_index() const				{return m_gridDataIndex;}

	protected:
		uint						m_gridDataIndex;//	index to grid-attached data.
};



////////////////////////////////////////////////////////////////////////////////////////////////
//	Vertex
///	Base-class for all vertex-types
/**
 * Vertices are required in any grid.
 * They are the geometric objects of lowest dimension.
 * All other geometric objects of higher dimension reference vertices.
 *
 * \ingroup lib_grid_grid_objects
 */
class UG_API Vertex : public GridObject
{
	friend class Grid;
	public:
		typedef Vertex grid_base_object;

	//	lower dimensional Base Object
		typedef void lower_dim_base_object;

	//	higher dimensional Base Object
		typedef Edge higher_dim_base_object;

	/**	The side type is obviously wrong. It should be void.
	 * However, void would cause problems with template instantiations.*/
		typedef Vertex side;
		typedef Edge sideof;

		static const bool HAS_SIDES = false;
		static const bool CAN_BE_SIDE = true;

	/// reference dimension
		static const int dim = 0;

		static const int BASE_OBJECT_ID = VERTEX;

		static const size_t NUM_VERTICES = 1;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<Vertex*>(pObj) != NULL;}

		virtual ~Vertex()	{}

		inline uint num_sides() const	{return 0;}

		virtual int container_section() const	{return -1;}
		virtual int base_object_id() const		{return VERTEX;}
		virtual ReferenceObjectID reference_object_id() const	{return ROID_UNKNOWN;}

	///	returns a value that can be used for hashing.
	/**	this value is unique per grid.
	 *	It is only set correctly if the vertex has been created
	 *	by the grid or has been properly registered at it.*/
		inline uint32 get_hash_value() const	{return m_hashValue;}

	protected:
		uint32	m_hashValue;//	a unique value for each vertex in a grid.
};


///	Base class for all classes which consist of a group of vertices
class UG_API IVertexGroup
{
	public:
		typedef Vertex* const* ConstVertexArray;

		virtual ~IVertexGroup()	{};
		virtual Vertex* vertex(uint index) const = 0;
		virtual ConstVertexArray vertices() const = 0;
		virtual size_t num_vertices() const = 0;

	//	compatibility with std::vector for some template routines
	///	returns the number of vertices.
		inline size_t size() const	{return num_vertices();}
	///	returns the i-th vertex.
		inline Vertex* operator[](size_t index) const {return vertex(index);}
};


///	this class can be used if one wants to create a custom element from a set of vertices.
class UG_API CustomVertexGroup : public IVertexGroup
{
	public:
		CustomVertexGroup()							{}
		CustomVertexGroup(size_t numVrts)			{m_vrts.resize(numVrts);}
		virtual ~CustomVertexGroup()				{}

		virtual Vertex* vertex(uint index) const	{return m_vrts[index];}
		virtual ConstVertexArray vertices() const	{return &m_vrts.front();}
		virtual size_t num_vertices() const			{return m_vrts.size();}

		void set_num_vertices(size_t newSize)		{m_vrts.resize(newSize);}
		void resize(size_t newSize)					{m_vrts.resize(newSize);}
		void clear()								{m_vrts.clear();}
		void push_back(Vertex* vrt)					{m_vrts.push_back(vrt);}
		void set_vertex(size_t index, Vertex* vrt)	{m_vrts[index] = vrt;}

	private:
		std::vector<Vertex*>	m_vrts;
};



////////////////////////////////////////////////////////////////////////////////////////////////
//	EdgeVertices
///	holds the vertices of an Edge or an EdgeDescriptor.
class UG_API EdgeVertices : public IVertexGroup
{
	friend class Grid;
	public:
		virtual ~EdgeVertices()						{}
		virtual Vertex* vertex(uint index) const	{return m_vertices[index];}
		virtual ConstVertexArray vertices() const	{return m_vertices;}
		virtual size_t num_vertices() const			{return 2;}

	//	compatibility with std::vector for some template routines
	///	returns the number of vertices.
		inline size_t size() const					{return 2;}
	///	returns the i-th vertex.
		Vertex* operator[](size_t index) const 		{return m_vertices[index];}

	protected:
		inline void assign_edge_vertices(const EdgeVertices& ev)
		{
			m_vertices[0] = ev.m_vertices[0];
			m_vertices[1] = ev.m_vertices[1];
		}

	protected:
		Vertex*	m_vertices[2];
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	Edge
///	Base-class for edges
/**
 * Edge is the base class of all 1-dimensional geometric objects.
 * Edges connect two vertices.
 *
 * \ingroup lib_grid_grid_objects
 */
class UG_API Edge : public GridObject, public EdgeVertices
{
	friend class Grid;
	public:
		typedef Edge grid_base_object;

	//	lower dimensional Base Object
		typedef Vertex lower_dim_base_object;
	//	higher dimensional Base Object
		typedef Face higher_dim_base_object;

		typedef Vertex side;
		typedef Face sideof;

		static const bool HAS_SIDES = true;
		static const bool CAN_BE_SIDE = true;

	/// reference dimension
		static const int dim = 1;

		static const int BASE_OBJECT_ID = EDGE;

		static const size_t NUM_VERTICES = 2;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<Edge*>(pObj) != NULL;}

		virtual ~Edge()	{}

		virtual int container_section() const	{return -1;}
		virtual int base_object_id() const		{return EDGE;}
		virtual ReferenceObjectID reference_object_id() const	{return ROID_UNKNOWN;}

		inline uint num_sides() const	{return 2;}

	///	retrieves the vertex on the opposing side to the specified one.
	/**	If the specified vertex is not part of the edge, false is returned.
	 * If it is, then vrtOut is filled with the opposing vertex and true is returned.*/
		bool get_opposing_side(Vertex* v, Vertex** vrtOut);

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
		virtual bool refine(std::vector<Edge*>& vNewEdgesOut,
											Vertex* newVertex,
											Vertex** pSubstituteVrts = NULL)	{return false;}

	protected:
		inline void set_vertex(uint index, Vertex* pVrt)	{m_vertices[index] = pVrt;}
};


////////////////////////////////////////////////////////////////////////////////////////////////
//	EdgeDescriptor
///	Can be used to store information about an edge and to construct an edge.
class UG_API EdgeDescriptor : public EdgeVertices
{
	public:
		EdgeDescriptor();
		EdgeDescriptor(const EdgeDescriptor& ed);
		EdgeDescriptor(Vertex* vrt1, Vertex* vrt2);

		EdgeDescriptor& operator = (const EdgeDescriptor& ed);

		inline void set_vertex(uint index, Vertex* vrt)	{m_vertices[index] = vrt;}
		inline void set_vertices(Vertex* vrt1, Vertex* vrt2)
			{
				m_vertices[0] = vrt1;
				m_vertices[1] = vrt2;
			}
};



class UG_API FaceVertices : public IVertexGroup
{
	public:
		virtual ~FaceVertices()							{}
		virtual Vertex* vertex(uint index) const		{UG_ASSERT(0, "SHOULDN'T BE CALLED"); return NULL;}
		virtual ConstVertexArray vertices() const		{UG_ASSERT(0, "SHOULDN'T BE CALLED"); return NULL;}
		virtual size_t num_vertices() const				{UG_ASSERT(0, "SHOULDN'T BE CALLED"); return 0;}

	//	compatibility with std::vector for some template routines
	///	returns the number of vertices.
		inline size_t size() const	{return num_vertices();}
	///	returns the i-th vertex.
		inline Vertex* operator[](size_t index) const {return vertex(index);}
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
 * \ingroup lib_grid_grid_objects
 */
class UG_API Face : public GridObject, public FaceVertices
{
	friend class Grid;
	public:
		typedef Face grid_base_object;

	//	lower dimensional Base Object
		typedef Edge lower_dim_base_object;
	//	higher dimensional Base Object
		typedef Volume higher_dim_base_object;

		typedef Edge side;
		typedef Volume sideof;

		static const bool HAS_SIDES = true;
		static const bool CAN_BE_SIDE = true;

	/// reference dimension
		static const int dim = 2;
		static const int BASE_OBJECT_ID = FACE;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<Face*>(pObj) != NULL;}

		virtual ~Face()	{}

	///	returns the i-th edge of the face.
	/**	This default implementation is reimplemented by derived classes for optimal speed.*/
		virtual EdgeDescriptor edge_desc(int index) const
			{return EdgeDescriptor(vertex(index), vertex((index+1) % size()));}

	///	returns the i-th edge of the face.
	/**	This default implementation is reimplemented by derived classes for optimal speed.*/
		virtual void edge_desc(int index, EdgeDescriptor& edOut) const
			{edOut.set_vertices(vertex(index), vertex((index+1) % size()));}

		inline uint num_edges() const	{return num_vertices();}
		inline uint num_sides() const	{return num_edges();}

		virtual int container_section() const	{return -1;}
		virtual int base_object_id() const		{return FACE;}
		virtual ReferenceObjectID reference_object_id() const	{return ROID_UNKNOWN;}

	/**	A default implementation is featured to allow empty instances of
	 *	this class. This is required to allow the use of this class
	 *	for compile-time method selection by dummy-parameters.
	 *	It is cruical that derived classes overload this method.*/
		virtual Edge* create_edge(int index)	{return NULL;}	///< create the edge with index i and return it.


	///	retrieves the edge-descriptor for the opposing side to the specified one.
	/**	If no opposing side exists false is returned. If an opposing side exists,
	 * the method returns true and fills the specified descriptor.*/
		virtual bool get_opposing_side(EdgeVertices* e, EdgeDescriptor& edOut) const	{return false;}

	///	returns an index to the obbject, which lies on the opposite side of the face to the given vertex.
	/**	The method returs a pair <GridBaseObjectId, int>, where GridBaseObjectId
	 * is either VERTEX or EDGE and where the second entry specifies the local index
	 * of the object in the given face.*/
		virtual std::pair<GridBaseObjectId, int>
		get_opposing_object(Vertex* vrt) const	{UG_THROW("Base method not implemented.");}

	///	returns the local index of the specified edge.
	/**	If the edge is not part of the face, then -1 is returned.*/
		int get_local_side_index(EdgeVertices* e) const;

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
							Vertex** newFaceVertexOut,
							Vertex** newEdgeVertices,
							Vertex* newFaceVertex = NULL,
							Vertex** pSubstituteVertices = NULL)	{return false;}

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
								int edgeIndex, Vertex* newVertex,
								Vertex** pSubstituteVertices = NULL)	{return false;}

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
								std::vector<Vertex*>& vNewEdgeVertices,
								Vertex** pSubstituteVertices = NULL)	{return false;}

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
							Vertex* newVertex,
							std::vector<Face*>& vNewFacesOut,
							Vertex** pSubstituteVertices = NULL)	{};
// END Depreciated

	/**	creates the faces that result from the collapsing of the edge with index 'splitEdgeIndex'.*/
		//virtual void create_faces_by_edge_collapse(int collapseEdgeIndex,
		//						Vertex* newVertex,
		//						std::vector<Face*>& vNewFacesOut) = 0;

	protected:
		virtual void set_vertex(uint index, Vertex* pVrt)	{UG_ASSERT(0, "SHOULDN'T BE CALLED");}
};


///	constant that defines the maximal number of vertices per face.
/**	This constant is mainly used by FaceDescriptor.
 * If required, one should be able to increase it without any problems.*/
const int MAX_FACE_VERTICES = 4;

////////////////////////////////////////////////////////////////////////////////////////////////
//	FaceDescriptor
///	Can be queried for the edges and vertices of a face.
class UG_API FaceDescriptor : public FaceVertices
{
	public:
		FaceDescriptor();
		FaceDescriptor(uint numVertices);
		FaceDescriptor(const FaceDescriptor& fd);
		FaceDescriptor(Vertex* v0, Vertex* v1, Vertex* v2) :
			m_numVertices(3)
			{m_vertices[0] = v0; m_vertices[1] = v1; m_vertices[2] = v2;}

		FaceDescriptor(Vertex* v0, Vertex* v1, Vertex* v2, Vertex* v3) :
			m_numVertices(4)
			{m_vertices[0] = v0; m_vertices[1] = v1; m_vertices[2] = v2; m_vertices[3] = v3;}

		virtual ~FaceDescriptor()					{}

		FaceDescriptor& operator = (const FaceDescriptor& fd);

		virtual Vertex* vertex(uint index) const	{return m_vertices[index];}
		virtual ConstVertexArray vertices() const		{return m_vertices;}
		virtual size_t num_vertices() const				{return m_numVertices;}

		inline void set_num_vertices(uint numVertices)	{m_numVertices = numVertices;}
		inline void set_vertex(uint index, Vertex* vrt)
			{m_vertices[index] = vrt;}

	protected:
		Vertex*	m_vertices[MAX_FACE_VERTICES];
		uint		m_numVertices;
};


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
//	VOLUMES

////////////////////////////////////////////////////////////////////////////////////////////////
//	VolumeVertices
///	holds the vertices of a Volume or a VolumeDescriptor
class UG_API VolumeVertices : public IVertexGroup
{
	public:
		virtual ~VolumeVertices()						{}

		virtual Vertex* vertex(uint index) const	{UG_ASSERT(0, "SHOULDN'T BE CALLED"); return NULL;}
		virtual ConstVertexArray vertices() const		{UG_ASSERT(0, "SHOULDN'T BE CALLED"); return NULL;}
		virtual size_t num_vertices() const				{UG_ASSERT(0, "SHOULDN'T BE CALLED"); return 0;}

	//	compatibility with std::vector for some template routines
	///	returns the number of vertices.
		inline size_t size() const	{return num_vertices();}
	///	returns the i-th vertex.
		inline Vertex* operator[](size_t index) const {return vertex(index);}
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
 * It is cruical that derived classes overload those methods.
 *
 * \ingroup lib_grid_grid_objects
 */
class UG_API Volume : public GridObject, public VolumeVertices
{
	friend class Grid;
	public:
		typedef Volume grid_base_object;

	//	lower dimensional Base Object
		typedef Face lower_dim_base_object;
	//	higher dimensional Base Object
		typedef void higher_dim_base_object;

		typedef Face side;
	//	this is just by convention. Use can_be_side() to stop recursions.
		typedef Volume sideof;

		static const bool HAS_SIDES = true;
		static const bool CAN_BE_SIDE = false;

	/// reference dimension
		static const int dim = 3;
		static const int BASE_OBJECT_ID = VOLUME;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<Volume*>(pObj) != NULL;}

		virtual ~Volume()	{}

		virtual EdgeDescriptor edge_desc(int index) const				{return EdgeDescriptor(NULL, NULL);}
		virtual void edge_desc(int index, EdgeDescriptor& edOut) const	{edOut = EdgeDescriptor(NULL, NULL);}
		virtual uint num_edges() const									{return 0;}

		virtual FaceDescriptor face_desc(int index) const				{return FaceDescriptor(0);}
		virtual void face_desc(int index, FaceDescriptor& fdOut) const	{fdOut = FaceDescriptor(0);}
		virtual uint num_faces() const									{return 0;}
		inline uint num_sides() const									{return num_faces();}

		virtual Edge* create_edge(int index)	{return NULL;}	///< create the edge with index i and return it.
		virtual Face* create_face(int index)		{return NULL;}	///< create the face with index i and return it.
		
	///	returns the local indices of an edge of the volume.
	/**	Default implementation throws an instance of int.
	 *	This should be changed by making the method pure virtual.*/
		virtual void get_local_vertex_indices_of_edge(size_t& ind1Out,
													  size_t& ind2Out,
													  size_t edgeInd) const
		{
		//	("Missing implementation of get_local_vertex_indices_of_face.")
			UG_THROW("Base method not implemented.");
		};
		
	///	returns the local indices of a face of the volume.
	/**	Default implementation throws an instance of int.
	 *	This should be changed by making the method pure virtual.*/
		virtual void get_local_vertex_indices_of_face(std::vector<size_t>& indsOut,
													  size_t side) const
		{
		//	("Missing implementation of get_local_vertex_indices_of_face.")
			UG_THROW("Base method not implemented.");
		};

	///	retrieves the face-descriptor for the opposing side to the specified one.
	/**	If no opposing side exists false is returned. If an opposing side exists,
	 * the method returns true and fills the specified descriptor.*/
		virtual bool get_opposing_side(FaceVertices* f, FaceDescriptor& fdOut) const	{return false;}

	///	returns an index to the obbject, which lies on the opposite side of the volume to the given vertex.
	/**	The method returs a pair <GridBaseObjectId, int>, where GridBaseObjectId
	 * is either VERTEX, EDGE, or FACE and where the second entry specifies the local index
	 * of the object in the given volume.*/
		virtual std::pair<GridBaseObjectId, int>
		get_opposing_object(Vertex* vrt) const	{UG_THROW("Base method not implemented.");}

	///	returns the local index of the given face or -1, if the face is not part of the volume.
		int get_local_side_index(FaceVertices* f) const;

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
	 *   an instance of a standard vertex (a concrete type is required. Vertex will not do).
	 *   The prototypeVertex has not to be registered at any grid - it may be a temporary
	 *   instance.
	 * - If you specify pvSubstituteVertices, the created volumes will reference the vertices
	 *   in pvSubstituteVertices. Note that pvSubstituteVertices has to contain exactly as
	 *   many vertices as the refined Face. Vertices with the same index correlate.
	 * - You may optionally specify an array containing the corner positions of the element.
	 * 	 If specified, those corners are used during regular tetrahedron refinement
	 * 	 only (only if all edges are refined). They are used to choose the best
	 * 	 refinement rule for interior child tetrahedrons, thus improving interior
	 * 	 angles.
	 */
		virtual bool refine(std::vector<Volume*>& vNewVolumesOut,
							Vertex** ppNewVertexOut,
							Vertex** newEdgeVertices,
							Vertex** newFaceVertices,
							Vertex* newVolumeVertex,
							const Vertex& prototypeVertex,
							Vertex** pSubstituteVertices = NULL,
							vector3* corners = NULL)	{return false;}

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
								int edgeIndex, Vertex* newVertex,
								std::vector<Vertex*>* pvSubstituteVertices = NULL)	{return false;}

	/**
	 * Writes vertices to the volume-descriptor so that it defines a volume with
	 * flipped orientation. If you want to flip the orientation of a volume in a
	 * grid, please consider using the grids flip_orientation method.
	 *
	 * Please note: The default implementation returns the original volume and has
	 * to be reimplemented by derived classes.
	 */
	 	virtual void get_flipped_orientation(VolumeDescriptor& vdOut) const;
		
		virtual int container_section() const	{return -1;}
		virtual int base_object_id() const		{return VOLUME;}
		virtual ReferenceObjectID reference_object_id() const	{return ROID_UNKNOWN;}

	/**	creates the volumes that result from the splitting of the edge with index 'splitEdgeIndex'.*/
		//virtual void create_volumes_by_edge_split(int splitEdgeIndex,
		//						Vertex* newVertex,
		//						std::vector<Volume*>& vNewFacesOut) = 0;

	/**	creates the volumes that result from the collapsing of the edge with index 'splitEdgeIndex'.*/
		//virtual void create_Volumes_by_edge_collapse(int collapseEdgeIndex,
		//						Vertex* newVertex,
		//						std::vector<Volume*>& vNewFacesOut) = 0;
	protected:
		virtual void set_vertex(uint index, Vertex* pVrt)	{UG_ASSERT(0, "SHOULDN'T BE CALLED");}
};


///	constant that defines the maximal number of vertices per volume element.
/**	This constant is mainly used by VolumeDescriptor.*/
const int MAX_VOLUME_VERTICES = 8;

////////////////////////////////////////////////////////////////////////////////////////////////
//	VolumeDescriptor
///	Holds a set of vertices which represent the corners of a volume element
class UG_API VolumeDescriptor : public VolumeVertices
{
	public:
		VolumeDescriptor();
		VolumeDescriptor(uint numVertices);
		VolumeDescriptor(const VolumeDescriptor& vd);

		virtual ~VolumeDescriptor()										{}

		VolumeDescriptor& operator = (const VolumeDescriptor& vv);
		VolumeDescriptor& operator = (const VolumeVertices& vv);

		virtual Vertex* vertex(uint index) const	{return m_vertices[index];}
		virtual ConstVertexArray vertices() const		{return m_vertices;}
		virtual size_t num_vertices() const				{return m_numVertices;}

		inline void set_num_vertices(uint numVertices)		{m_numVertices = numVertices;}
		inline void set_vertex(uint index, Vertex* vrt)	{m_vertices[index] = vrt;}

	protected:
		Vertex*	m_vertices[MAX_VOLUME_VERTICES];
		uint		m_numVertices;
};



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/** Defines the geometric base object type for each dimension.
 * \{ */
template <int dim> struct GeomObjBaseTypeByDim;

template <> struct GeomObjBaseTypeByDim<0>{
	typedef Vertex base_obj_type;
};

template <> struct GeomObjBaseTypeByDim<1>{
	typedef Edge base_obj_type;
};

template <> struct GeomObjBaseTypeByDim<2>{
	typedef Face base_obj_type;
};

template <> struct GeomObjBaseTypeByDim<3>{
	typedef Volume base_obj_type;
};

/** \} */


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/**	template helpers that return the geometric base object type
 *	given a pointer to a derived class of Vertex, Edge, Face or Volume.
 *
 *	e.g. PtrTypeToGeomObjBaseType<RegularVertex*>::base_type = Vertex.
 * \{
 */
template <class TGeomObjPtrType>
struct PtrTypeToGeomObjBaseType
{typedef void base_type;};

template <>
struct PtrTypeToGeomObjBaseType<Vertex*>
{typedef Vertex base_type;};

template <>
struct PtrTypeToGeomObjBaseType<Edge*>
{typedef Edge base_type;};

template <>
struct PtrTypeToGeomObjBaseType<Face*>
{typedef Face base_type;};

template <>
struct PtrTypeToGeomObjBaseType<Volume*>
{typedef Volume base_type;};
/** \} */


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	hash-funtions for vertices
///	returns the hash-value of the vertex.
template <>
size_t hash_key<PVertex>(const PVertex& key);

////////////////////////////////////////////////////////////////////////
//	hash-funtions for edges
///	the hash-key is a function of vertex-hash-values.
/**
 * The hash value depends on the associated vertices.
 * If an Edge (or EdgeDescriptor) has the same vertices
 * as another Edge (or EdgeDescriptor), the hash-keys
 * are the same.
 */
template <>
size_t hash_key<PEdgeVertices>(const PEdgeVertices& key);

///	the hash-key is a function of vertex-hash-values.
/** \sa hash_key<PEdgeVertices>*/
template <>
size_t hash_key<PEdge>(const PEdge& key);

///	the hash-key is a function of vertex-hash-values.
/** \sa hash_key<PEdgeVertices>*/
template <>
size_t hash_key<PEdgeDescriptor>(const PEdgeDescriptor& key);

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
size_t hash_key<PFaceVertices>(const PFaceVertices& key);

///	the hash-key is a function of vertex-hash-values.
/**\sa hash_key<PFaceVertices>*/
template <>
size_t hash_key<PFace>(const PFace& key);

///	the hash-key is a function of vertex-hash-values.
/**\sa hash_key<PFaceVertices>*/
template <>
size_t hash_key<PFaceDescriptor>(const PFaceDescriptor& key);

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
size_t hash_key<PVolumeVertices>(const PVolumeVertices& key);

///	the hash-key is a function of vertex-hash-values.
/**\sa hash_key<PVolumeVertices>*/
template <>
size_t hash_key<PVolume>(const PVolume& key);

///	the hash-key is a function of vertex-hash-values.
/**\sa hash_key<PVolumeVertices>*/
template <>
size_t hash_key<PVolumeDescriptor>(const PVolumeDescriptor& key);

}//	end of namespace

#endif
