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

#ifndef __H__LIB_GRID__GEOMETRIC_BASE_OBJECTS__
#define __H__LIB_GRID__GEOMETRIC_BASE_OBJECTS__

//#include <list>
#include <cassert>
//#include <iostream>
#include <utility>
#include "common/types.h"
#include "common/common.h"
#include "common/assert.h"
#include "lib_grid/attachments/attachment_pipe.h"
#include "lib_grid/attachments/attached_list.h"
//#include "common/util/hash_function.h"
#include "common/math/ugmath_types.h"
//#include "common/util/pointer_const_array.h"

namespace ug {

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	CONSTANTS

////////////////////////////////////////////////////////////////////////
///	enumeration of the GridBaseObjects that make up a grid.
using GridBaseObjectId_t = uint;
enum GridBaseObjectId : uint
{
	VERTEX = 0,
	EDGE,
	FACE,
	VOLUME,
	NUM_GEOMETRIC_BASE_OBJECTS		//always last!!!
};

extern const char* GRID_BASE_OBJECT_SINGULAR_NAMES[];
extern const char* GRID_BASE_OBJECT_PLURAL_NAMES[];

////////////////////////////////////////////////////////////////////////
//	Reference-Object IDs
///	these ids are used to identify the shape of a geometric object.

enum ReferenceObjectID : int
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
using ReferenceObjectID_t = ReferenceObjectID;

inline
ReferenceObjectID operator ++ (ReferenceObjectID& roid, int)
{
	int tmp = roid;
	return roid = static_cast<ReferenceObjectID>(++tmp);
};

////////////////////////////////////////////////////////////////////////
inline
std::ostream& operator << (std::ostream& outStream, ReferenceObjectID type)
{
	switch(type)
	{
		case ReferenceObjectID::ROID_UNKNOWN: outStream << "(invalid)"; break;
		case ReferenceObjectID::ROID_VERTEX: outStream << "Vertex"; break;
		case ReferenceObjectID::ROID_EDGE: outStream << "Edge"; break;
		case ReferenceObjectID::ROID_TRIANGLE: outStream << "Triangle"; break;
		case ReferenceObjectID::ROID_QUADRILATERAL: outStream << "Quadrilateral"; break;
		case ReferenceObjectID::ROID_TETRAHEDRON: outStream << "Tetrahedron"; break;
		case ReferenceObjectID::ROID_OCTAHEDRON: outStream << "Octahedron"; break;
		case ReferenceObjectID::ROID_HEXAHEDRON: outStream << "Hexahedron"; break;
		case ReferenceObjectID::ROID_PRISM: outStream << "Prism"; break;
		case ReferenceObjectID::ROID_PYRAMID: outStream << "Pyramid"; break;
		default: throw(UGError("Unknown ReferenceObjectID in operator <<"));
	}
	return outStream;
};


////////////////////////////////////////////////////////////////////////
//	PREDECLARATIONS
class Grid;
template <typename TElem>	class ElementStorage;

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
class UG_API GridObject
{
	friend class Grid;
	friend class attachment_traits<Vertex*, ElementStorage<Vertex> >;
	friend class attachment_traits<Edge*, ElementStorage<Edge> >;
	friend class attachment_traits<Face*, ElementStorage<Face> >;
	friend class attachment_traits<Volume*, ElementStorage<Volume> >;

	public:
		virtual ~GridObject() = default;

	///	create an instance of the derived type
	/**	Make sure to overload this method in derivates of this class!*/
		[[nodiscard]] virtual GridObject* create_empty_instance() const {return nullptr;}

		[[nodiscard]] virtual int container_section() const = 0;
		[[nodiscard]] virtual GridBaseObjectId_t base_object_id() const = 0;
	/**
	 * A reference object represents a class of geometric objects.
	 * Tetrahedrons, Triangles etc are such classes.
	 * Reference ids should be defined in the file in which concrete geometric objects are defined.
	 */
		[[nodiscard]] virtual ReferenceObjectID_t reference_object_id() const = 0;///	returns the id of the reference-object.

	///	returns true if the object constrains other objects.
	/**	This is normally only the case for special constraining objects.
	 * The default implementation returns false.*/
		[[nodiscard]] virtual bool is_constraining() const {return false;}

	///	returns true if the object is constrained by other objects.
	/**	This is normally only the case for special constrained objects.
	 * The default implementation returns false.*/
		[[nodiscard]] virtual bool is_constrained() const {return false;}

	///	removes a constraint link to the grid object.
	/**	This method is e.g. called on the constraining edges of a constrained
	 * vertex as soon as the constrained vertex is deleted.
	 * \{ */
		virtual void remove_constraint_link(const Vertex* vrt)		{}
		virtual void remove_constraint_link(const Edge* e)			{}
		virtual void remove_constraint_link(const Face* f)			{}
	/** \} */

	///	Returns the grid attachment data index of a geometric object.
	/** Beware that this index is for internal use in the grid management and can
	 * be changed by some operations (e.g. by Grid::defragment). But this function
	 * can be used for debugging, when one wants to identify an element at several
	 * places in the code.
	 */
		[[nodiscard]] inline uint grid_data_index() const {return m_gridDataIndex;}

	protected:
	///	ATTENTION: Use this method with extreme care!
	/**	This method is for internal use only and should almost never be called
	 * by a user of lib_grid. The method sets the attachment data index and is
	 * mainly used by attachment-traits classes.*/
		inline void set_grid_data_index(uint index) {m_gridDataIndex = index;}

	protected:
		uint m_gridDataIndex;//	index to grid-attached data.
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
		using grid_base_object = Vertex;

	//	lower dimensional Base Object
		using lower_dim_base_object = void;

	//	higher dimensional Base Object
		using higher_dim_base_object = Edge;

	/**	The side type is obviously wrong. It should be void.
	 * However, void would cause problems with template instantiations.*/
		using side = Vertex;
		using sideof = Edge;

		static constexpr bool HAS_SIDES = false;
		static constexpr bool CAN_BE_SIDE = true;

	/// reference dimension
		static constexpr int dim = 0;

		static constexpr int BASE_OBJECT_ID = GridBaseObjectId::VERTEX;

		static constexpr size_t NUM_VERTICES = 1;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<Vertex*>(pObj) != nullptr;}

		~Vertex() override = default;

		[[nodiscard]] inline uint num_sides() const {return 0;}

		[[nodiscard]] int container_section() const override {return -1;}
		[[nodiscard]] GridBaseObjectId_t base_object_id() const override {return GridBaseObjectId::VERTEX;}
		[[nodiscard]] ReferenceObjectID_t reference_object_id() const override {return ReferenceObjectID::ROID_UNKNOWN;}

	///	returns a value that can be used for hashing.
	/**	this value is unique per grid.
	 *	It is only set correctly if the vertex has been created
	 *	by the grid or has been properly registered at it.*/
		[[nodiscard]] inline uint32 get_hash_value() const	{return m_hashValue;}

	protected:
		uint32	m_hashValue;//	a unique value for each vertex in a grid.
};


///	This descriptor is mainly useful to avoid compilation errors in templated code
/** \note	Most methods here are only provided to avoid issues with other descriptors
 *			in templated code.
 *	\note	The only valid index is '0' for any method taking indices. If a higher
 *			index is passed, the index is ignored and '0' is used instead.
 */
class UG_API VertexDescriptor {
public:
	using ConstVertexArray = Vertex* const*;

	VertexDescriptor() = default;
	VertexDescriptor(Vertex* v) : m_v(v)	{}
	VertexDescriptor(const VertexDescriptor& d) : m_v (d.m_v)	{}
	virtual ~VertexDescriptor() = default;

	VertexDescriptor& operator = (const VertexDescriptor& d)
	{m_v = d.m_v; return *this;}

	inline void set_vertex(Vertex* v)		{m_v = v;}
	inline void set_vertex(uint, Vertex* v)	{m_v = v;}

	[[nodiscard]] inline Vertex* vertex () const {return m_v;}
	[[nodiscard]] virtual Vertex* vertex(size_t) const {return m_v;}
	[[nodiscard]] Vertex* operator [] (size_t) const {return m_v;}

	[[nodiscard]] virtual ConstVertexArray vertices() const	{return &m_v;}
	[[nodiscard]] virtual size_t num_vertices() const {return 1;}

//	compatibility with std::vector for some template routines
///	returns the number of vertices.
	[[nodiscard]] inline size_t size() const    {return 1;}

private:
	Vertex* m_v;
};


///	Base class for all classes which consist of a group of vertices
class UG_API IVertexGroup
{
	public:
		using ConstVertexArray = Vertex* const*;

		virtual ~IVertexGroup()	= default;
		virtual Vertex* vertex(size_t index) const = 0;
		virtual ConstVertexArray vertices() const = 0;
		virtual size_t num_vertices() const = 0;

	//	compatibility with std::vector for some template routines
	///	returns the number of vertices.
		[[nodiscard]] inline size_t size() const {return num_vertices();}
	///	returns the i-th vertex.
		inline Vertex* operator [] (size_t index) const {return vertex(index);}
};


///	this class can be used if one wants to create a custom element from a set of vertices.
class UG_API CustomVertexGroup : public IVertexGroup
{
	public:
		CustomVertexGroup() = default;
		explicit CustomVertexGroup(size_t numVrts) {m_vrts.resize(numVrts);}
		~CustomVertexGroup() override = default;

		[[nodiscard]] Vertex* vertex(size_t index) const override {return m_vrts[index];}
		[[nodiscard]] ConstVertexArray vertices() const override {return &m_vrts.front();}
		[[nodiscard]] size_t num_vertices() const override {return m_vrts.size();}

		void set_num_vertices(size_t newSize) {m_vrts.resize(newSize);}
		void resize(size_t newSize) {m_vrts.resize(newSize);}
		void clear() {m_vrts.clear();}
		void push_back(Vertex* vrt) {m_vrts.push_back(vrt);}
		void set_vertex(size_t index, Vertex* vrt) {m_vrts[index] = vrt;}

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
		~EdgeVertices() override = default;
		[[nodiscard]] Vertex* vertex(size_t index) const override {return m_vertices[index];}
		[[nodiscard]] ConstVertexArray vertices() const override {return m_vertices;}
		[[nodiscard]] size_t num_vertices() const override {return 2;}

	//	compatibility with std::vector for some template routines
	///	returns the number of vertices.
		[[nodiscard]] inline size_t size() const {return 2;}
	///	returns the i-th vertex.
		Vertex* operator [] (size_t index) const {return m_vertices[index];}

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
		using grid_base_object = Edge;

	//	lower dimensional Base Object
		using lower_dim_base_object = Vertex;
	//	higher dimensional Base Object
		using higher_dim_base_object = Face;

		using side = Vertex;
		using sideof = Face;

		static constexpr bool HAS_SIDES = true;
		static constexpr bool CAN_BE_SIDE = true;

	/// reference dimension
		static constexpr int dim = 1;

		static constexpr int BASE_OBJECT_ID = GridBaseObjectId::EDGE;

		static constexpr size_t NUM_VERTICES = 2;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<Edge*>(pObj) != nullptr;}

		~Edge() override = default;

		[[nodiscard]] int container_section() const override {return -1;}
		[[nodiscard]] GridBaseObjectId_t base_object_id() const override {return GridBaseObjectId::EDGE;}
		[[nodiscard]] ReferenceObjectID_t reference_object_id() const override {return ReferenceObjectID::ROID_UNKNOWN;}

		[[nodiscard]] inline uint num_sides() const {return 2;}

	///	retrieves the vertex on the opposing side to the specified one.
	/**	If the specified vertex is not part of the edge, false is returned.
	 * If it is, then vrtOut is filled with the opposing vertex and true is returned.*/
		bool get_opposing_side(Vertex* v, Vertex** vrtOut);

		bool get_opposing_side(Vertex* v, VertexDescriptor& vrtOut);

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
											Vertex** pSubstituteVrts = nullptr)	{return false;}

	protected:
		inline void set_vertex(uint index, Vertex* pVrt) {m_vertices[index] = pVrt;}
};


////////////////////////////////////////////////////////////////////////////////////////////////
//	EdgeDescriptor
///	Can be used to store information about an edge and to construct an edge.
class UG_API EdgeDescriptor : public EdgeVertices
{
	public:
		EdgeDescriptor() = default;
		EdgeDescriptor(const EdgeDescriptor& ed);
		EdgeDescriptor(Vertex* vrt1, Vertex* vrt2);
		~EdgeDescriptor() override = default;

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
		~FaceVertices() override = default;
		[[nodiscard]] Vertex* vertex(size_t index) const override {UG_ASSERT(0, "SHOULDN'T BE CALLED"); return nullptr;}
		[[nodiscard]] ConstVertexArray vertices() const override {UG_ASSERT(0, "SHOULDN'T BE CALLED"); return nullptr;}
		[[nodiscard]] size_t num_vertices() const override {UG_ASSERT(0, "SHOULDN'T BE CALLED"); return 0;}

	//	compatibility with std::vector for some template routines
	///	returns the number of vertices.
		[[nodiscard]] inline size_t size() const {return num_vertices();}
	///	returns the i-th vertex.
		inline Vertex* operator [] (size_t index) const {return vertex(index);}
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
		using grid_base_object = Face;

	//	lower dimensional Base Object
		using lower_dim_base_object = Edge;
	//	higher dimensional Base Object
		using higher_dim_base_object = Volume;

		using side = Edge;
		using sideof = Volume;

		static constexpr bool HAS_SIDES = true;
		static constexpr bool CAN_BE_SIDE = true;

	/// reference dimension
		static constexpr int dim = 2;
		static constexpr int BASE_OBJECT_ID = GridBaseObjectId::FACE;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<Face*>(pObj) != nullptr;}

		~Face() override = default;

	///	returns the i-th edge of the face.
	/**	This default implementation is reimplemented by derived classes for optimal speed.*/
		[[nodiscard]] virtual EdgeDescriptor edge_desc(int index) const
			{return EdgeDescriptor(vertex(index), vertex((index+1) % size()));}

	///	returns the i-th edge of the face.
	/**	This default implementation is reimplemented by derived classes for optimal speed.*/
		virtual void edge_desc(int index, EdgeDescriptor& edOut) const
			{edOut.set_vertices(vertex(index), vertex((index+1) % size()));}

		[[nodiscard]] inline uint num_edges() const	{return (uint)num_vertices();}
		[[nodiscard]] inline uint num_sides() const	{return num_edges();}

		[[nodiscard]] int container_section() const override {return -1;}
		[[nodiscard]] GridBaseObjectId_t base_object_id() const override {return GridBaseObjectId::FACE;}
		[[nodiscard]] ReferenceObjectID_t reference_object_id() const override {return ReferenceObjectID::ROID_UNKNOWN;}

	/**	A default implementation is featured to allow empty instances of
	 *	this class. This is required to allow the use of this class
	 *	for compile-time method selection by dummy-parameters.
	 *	It is cruical that derived classes overload this method.*/
		virtual Edge* create_edge(int index) {return nullptr;}	///< create the edge with index i and return it.


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
	 * are inserted on the edge that corresponds to their index. Use nullptr to indicate
	 * that no vertex shall be inserted on the associated edge. newEdgeVertices has to point
	 * to an array that holds as many vertices as there are edges in the face.
	 * - If the method has to create a new inner vertex, it will be returned in newFaceVertexOut.
	 * - If you specify pvSubstituteVertices, the created faces will reference the vertices in
	 * pvSubstituteVertices. Note that pvSubstituteVertices has to contain exactly as many
	 * vertices as the refined Face. Vertices with the same index correlate.
	 * - You may optionally specify a 'snapPointIndex'. Default is -1. If an index >= 0
	 *   is specifed, then new inner edges will always be created between new edge vertices
	 *   and the vertex specified through the snap-point-index. Note that a snap-point
	 *   must not be a corner of a refined edge.
	 */
		virtual bool refine(std::vector<Face*>& vNewFacesOut,
							Vertex** newFaceVertexOut,
							Vertex** newEdgeVertices,
							Vertex* newFaceVertex = nullptr,
							Vertex** pSubstituteVertices = nullptr,
							int snapPointIndex = -1)	{return false;}


	///	returns true if the specified edgeMarks would lead to a regular refinement
	/**	A regular refinement leads to new elements which are all similar to the
	 * original element. I.e. which are of the same type and which have similar
	 * angles.
	 *
	 * \note	this method does not perform refinement. Use 'refine' with the
	 *			specified edges instead.
	 *
	 * \param edgeMarks	If the i-th edge shall be refined, the expression
	 *					'edgeMarks & (1<<i) != 0' has to be true. You can
	 *					specify multiple refine-edges using or-combinations:
	 *					'edgeMarks = (1<<i) | (1<<j)' would indicate that the
	 *					i-th and the j-th edge shall be refined.*/
		[[nodiscard]] virtual bool is_regular_ref_rule(int edgeMarks) const	{return false;}

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
								Vertex** pSubstituteVertices = nullptr)	{return false;}

	/**
	 * The collapse_edgea method creates new geometric objects by collapsing the specified edges
	 * simultaneously. This method makes sense only for faces with more than 4 edges.
	 * The user that calls this function is responsible to either register the new
	 * faces with a grid (the grid from which the vertices are), or to take responsibility
	 * for deletion of the acquired memory (delete each element in vNewFacesOut).
	 * - for each entry in vNewEdgeVertices which is not nullptr, the edge with the same index will
	 *    be collapsed and replaced by the specified vertex.
	 *    The size of vNewEdgeVertices has thus to be between 0 and this->num_edges().
	 * - If you specify pvSubstituteVertices, the created faces will reference the vertices in
	 * pvSubstituteVertices. Note that pvSubstituteVertices has to contain exactly as many
	 * vertices as the face in which you collapse the edge. Vertices with the same index correlate.
	 */
		virtual bool collapse_edges(std::vector<Face*>& vNewFacesOut,
								std::vector<Vertex*>& vNewEdgeVertices,
								Vertex** pSubstituteVertices = nullptr)	{return false;}

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
							Vertex** pSubstituteVertices = nullptr)	{};
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
		explicit FaceDescriptor(uint numVertices);
		FaceDescriptor(const FaceDescriptor& fd);
		FaceDescriptor(Vertex* v0, Vertex* v1, Vertex* v2) :
			m_numVertices(3)
			{m_vertices[0] = v0; m_vertices[1] = v1; m_vertices[2] = v2;}

		FaceDescriptor(Vertex* v0, Vertex* v1, Vertex* v2, Vertex* v3) :
			m_numVertices(4)
			{m_vertices[0] = v0; m_vertices[1] = v1; m_vertices[2] = v2; m_vertices[3] = v3;}

		~FaceDescriptor() override = default;

		FaceDescriptor& operator = (const FaceDescriptor& fd);

		[[nodiscard]] Vertex* vertex(size_t index) const override {return m_vertices[index];}
		[[nodiscard]] ConstVertexArray vertices() const override {return m_vertices;}
		[[nodiscard]] size_t num_vertices() const override {return m_numVertices;}

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
		~VolumeVertices() override = default;

		[[nodiscard]] Vertex* vertex(size_t index) const override {UG_ASSERT(0, "SHOULDN'T BE CALLED"); return nullptr;}
		[[nodiscard]] ConstVertexArray vertices() const override {UG_ASSERT(0, "SHOULDN'T BE CALLED"); return nullptr;}
		[[nodiscard]] size_t num_vertices() const override {UG_ASSERT(0, "SHOULDN'T BE CALLED"); return 0;}

	//	compatibility with std::vector for some template routines
	///	returns the number of vertices.
		[[nodiscard]] inline size_t size() const {return num_vertices();}
	///	returns the i-th vertex.
		inline Vertex* operator [] (size_t index) const {return vertex(index);}
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
		using grid_base_object = Volume;

	//	lower dimensional Base Object
		using lower_dim_base_object = Face;
	//	higher dimensional Base Object
		using higher_dim_base_object = void;

		using side = Face;
	//	this is just by convention. Use can_be_side() to stop recursions.
		using sideof = Volume;

		static constexpr bool HAS_SIDES = true;
		static constexpr bool CAN_BE_SIDE = false;

	/// reference dimension
		static constexpr int dim = 3;
		static constexpr int BASE_OBJECT_ID = GridBaseObjectId::VOLUME;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<Volume*>(pObj) != nullptr;}

		~Volume() override = default;

		[[nodiscard]] virtual EdgeDescriptor edge_desc(int index) const {return EdgeDescriptor(nullptr, nullptr);}
		virtual void edge_desc(int index, EdgeDescriptor& edOut) const {edOut = EdgeDescriptor(nullptr, nullptr);}
		[[nodiscard]] virtual uint num_edges() const {return 0;}

		[[nodiscard]] virtual FaceDescriptor face_desc(int index) const {return FaceDescriptor(0);}
		virtual void face_desc(int index, FaceDescriptor& fdOut) const {fdOut = FaceDescriptor(0);}
		[[nodiscard]] virtual uint num_faces() const {return 0;}
		[[nodiscard]] inline uint num_sides() const {return num_faces();}

		virtual Edge* create_edge(int index) {return nullptr;}	///< create the edge with index i and return it.
		virtual Face* create_face(int index) {return nullptr;}	///< create the face with index i and return it.
		
	///	returns the local indices of an edge of the volume.
	/**	Default implementation throws a UGError*/
		virtual void get_vertex_indices_of_edge(size_t& ind1Out,
												size_t& ind2Out,
												size_t edgeInd) const
		{UG_THROW("Base method not implemented.");};
		
	///	returns the local indices of a face of the volume.
	/**	Default implementation throws a UGError*/
		virtual void get_vertex_indices_of_face(std::vector<size_t>& indsOut,
												size_t side) const
		{UG_THROW("Base method not implemented.");}

	///	returns the local index of the edge which connects the two vertex indices.
	/**	Default implementation throws a UGError*/
		virtual int get_edge_index_from_vertices(	const size_t vi0,
													const size_t vi1) const
		{UG_THROW("Base method not implemented.");};

	///	returns the local index of the j-th edge of the i-th face of the volume
	/**	Default implementation throws a UGError*/
		virtual int get_face_edge_index (	const size_t faceInd,
											const size_t faceEdgeInd) const
		{UG_THROW("Base method not implemented.");};

	///	retrieves the face-descriptor for the opposing side to the specified one.
	/**	If no opposing side exists false is returned. If an opposing side exists,
	 * the method returns true and fills the specified descriptor.*/
		virtual bool get_opposing_side(FaceVertices* f, FaceDescriptor& fdOut) const	{return false;}

	///	returns an index to the obbject, which lies on the opposite side of the volume to the given vertex.
	/**	The method returs a pair <GridBaseObjectId, int>, where GridBaseObjectId
	 * is either VERTEX, EDGE, or FACE and where the second entry specifies the local index
	 * of the object in the given volume.*/
		virtual std::pair<GridBaseObjectId, int>
		get_opposing_object(Vertex* vrt) const
		{UG_THROW("Base method not implemented.");}

	///	returns the local index of the given face or -1, if the face is not part of the volume.
		int get_local_side_index(FaceVertices* f) const;

	/**
	 * The refine method can be used to create new elements by inserting new vertices
	 * on the volume.
	 *
	 * New volumes will be returned in vNewFacesOut.
	 * If a new vertex has to be created from the prototypeVertex (this happens in
	 * more complicated situations) the pointer to the new vertex is returned in
	 * ppNewVertexOut. ppNewVertexOut contains nullptr if no new vertex has been created.
	 *
	 * The user that calls this function is responsible to either register the new
	 * volumes and the new vertex with a grid (the grid from which the referenced
	 * vertices are), or to take responsibility for deletion of the acquired memory
	 * (delete each element in vNewVolumesOut and ppNewVertexOut - if it is not nullptr).
	 *
	 * - Specify vertices that shall be inserted on edges with vNewEdgeVertices. Vertices
	 *   are inserted on the edge that corresponds to its index. Use nullptr to indicate
	 *   that no vertex shall be inserted on the associated edge.
	 * - Specify vertices that shall be inserted on faces with vNewFaceVertices. Vertices
	 *   are inserted on the face that corresponds to its index. Use nullptr to indicate
	 *   that no vertex shall be inserted on the associated face.
	 * - If newVolumeVertex is not nullptr, newVolumeVertex will be inserted in the center of
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
	 * - You may optionally specify 'isSnapPoint', a list of booleans of length 'num_vertices()'.
	 *   Each quadrilateral side may contain at most one snap-point. New edges on quadrilateral
	 *   faces will then connect the snap-point and the newly introduced edge-vertex.
	 *   Note that a snap-point must not be a corner of a refined edge.
	 */
		virtual bool refine(std::vector<Volume*>& vNewVolumesOut,
							Vertex** ppNewVertexOut,
							Vertex** newEdgeVertices,
							Vertex** newFaceVertices,
							Vertex* newVolumeVertex,
							const Vertex& prototypeVertex,
							Vertex** pSubstituteVertices = nullptr,
							vector3* corners = nullptr,
							bool* isSnapPoint = nullptr)	{return false;}

		
	///	returns true if the specified edgeMarks would lead to a regular refinement
	/**	A regular refinement leads to new elements which are all similar to the
	 * original element. I.e. which are of the same type and which have similar
	 * angles.
	 *
	 * \note	this method does not perform refinement. Use 'refine' with the
	 *			specified edges instead.
	 *
	 * \param edgeMarks	If the i-th edge shall be refined, the expression
	 *					'edgeMarks & (1<<i) != 0' has to be true. You can
	 *					specify multiple refine-edges using or-combinations:
	 *					'edgeMarks = (1<<i) | (1<<j)' would indicate that the
	 *					i-th and the j-th edge shall be refined.*/
		virtual bool is_regular_ref_rule(int edgeMarks) const {return false;}

	/**
	 * The collapse_edge method creates new geometric objects by collapsing the specified edge.
	 * The user that calls this function is responsible to either register the new
	 * volumes with a grid (the grid from which the vertices are), or to take responsibility
	 * for deletion of the acquired memory (delete each element in vNewFacesOut).
	 * - edgeIndex specifies the edge which shall be collapsed. edgeIndex has to lie
	 *   between 0 and this->num_edges().
	 * - Vertices adjacent to the collapsed edge will be replaced by newVertex.
	 * - If you specify pvSubstituteVertices, the created volumes will reference the vertices in
	 * pvSubstituteVertices. Note that pvSubstituteVertices has to contain exactly as many
	 * vertices as the volume in which you collapse the edge. Vertices with the same index correlate.
	 */
		virtual bool collapse_edge(std::vector<Volume*>& vNewVolumesOut,
								int edgeIndex, Vertex* newVertex,
								std::vector<Vertex*>* pvSubstituteVertices = nullptr)	{return false;}

	/**
	 * Writes vertices to the volume-descriptor so that it defines a volume with
	 * flipped orientation. If you want to flip the orientation of a volume in a
	 * grid, please consider using the grids flip_orientation method.
	 *
	 * Please note: The default implementation returns the original volume and has
	 * to be reimplemented by derived classes.
	 */
	 	virtual void get_flipped_orientation(VolumeDescriptor& vdOut) const;

		[[nodiscard]] int container_section() const override {return -1;}
		[[nodiscard]] GridBaseObjectId_t base_object_id() const override {return GridBaseObjectId::VOLUME;}
		[[nodiscard]] ReferenceObjectID_t reference_object_id() const override {return ReferenceObjectID::ROID_UNKNOWN;}

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
constexpr int MAX_VOLUME_VERTICES = 8;

////////////////////////////////////////////////////////////////////////////////////////////////
//	VolumeDescriptor
///	Holds a set of vertices which represent the corners of a volume element
class UG_API VolumeDescriptor : public VolumeVertices
{
	public:
		VolumeDescriptor() = default;
		explicit VolumeDescriptor(uint numVertices);
		VolumeDescriptor(const VolumeDescriptor& vd);

		~VolumeDescriptor() override = default;

		VolumeDescriptor& operator = (const VolumeDescriptor& vv);
		VolumeDescriptor& operator = (const VolumeVertices& vv);

		[[nodiscard]] Vertex* vertex(size_t index) const override {return m_vertices[index];}
		[[nodiscard]] ConstVertexArray vertices() const override {return m_vertices;}
		[[nodiscard]] size_t num_vertices() const override {return m_numVertices;}

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
	using base_obj_type = Vertex;
};

template <> struct GeomObjBaseTypeByDim<1>{
	using base_obj_type = Edge;
};

template <> struct GeomObjBaseTypeByDim<2>{
	using base_obj_type = Face;
};

template <> struct GeomObjBaseTypeByDim<3>{
	using base_obj_type = Volume;
};

/** \} */


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/**	template helpers that return the geometric base object type
 *	given a pointer to a derived class of Vertex, Edge, Face or Volume.
 *
 *	e.g. PtrToValueType<RegularVertex*>::base_type = Vertex.
 * \{
 */
template <typename TGeomObjPtrType>
struct PtrToValueType
{using base_type = void;};

template <>
struct PtrToValueType<Vertex*>
{using base_type = Vertex;};

template <>
struct PtrToValueType<Edge*>
{using base_type = Edge;};

template <>
struct PtrToValueType<Face*>
{using base_type = Face;};

template <>
struct PtrToValueType<Volume*>
{using base_type = Volume;};
/** \} */


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	hash-funtions for vertices
///	returns the hash-value of the vertex.
size_t hash_key(Vertex* key);

////////////////////////////////////////////////////////////////////////
//	hash-funtions for edges
///	the hash-key is a function of vertex-hash-values.
/**
 * The hash value depends on the associated vertices.
 * If an Edge (or EdgeDescriptor) has the same vertices
 * as another Edge (or EdgeDescriptor), the hash-keys
 * are the same.
 */
size_t hash_key(EdgeVertices* key);
size_t hash_key(const EdgeVertices* key);

///	the hash-key is a function of vertex-hash-values.
/** \sa hash_key<PEdgeVertices>*/
size_t hash_key(Edge* key);

///	the hash-key is a function of vertex-hash-values.
/** \sa hash_key<PEdgeVertices>*/
size_t hash_key(EdgeDescriptor* key);

////////////////////////////////////////////////////////////////////////
//	hash-funtions for faces
///	the hash-key is a function of vertex-hash-values.
/**
 * The hash value depends on the associated vertices.
 * If an Face (or FaceDescriptor) has the same vertices
 * as another Face (or FaceDescriptor), the hash-keys
 * are the same.
 */
size_t hash_key(FaceVertices* key);
size_t hash_key(const FaceVertices* key);

///	the hash-key is a function of vertex-hash-values.
/**\sa hash_key<PFaceVertices>*/
size_t hash_key(Face* key);

///	the hash-key is a function of vertex-hash-values.
/**\sa hash_key<PFaceVertices>*/
size_t hash_key(FaceDescriptor* key);

////////////////////////////////////////////////////////////////////////
//	hash-funtions for volumes
///	the hash-key is a function of vertex-hash-values.
/**
 * The hash value depends on the associated vertices.
 * If an Volume (or VolumeDescriptor) has the same vertices
 * as another Volume (or VolumeDescriptor), the hash-keys
 * are the same.
 */
size_t hash_key(VolumeVertices* key);
size_t hash_key(const VolumeVertices* key);

///	the hash-key is a function of vertex-hash-values.
/**\sa hash_key<PVolumeVertices>*/
size_t hash_key(Volume* key);

///	the hash-key is a function of vertex-hash-values.
/**\sa hash_key<PVolumeVertices>*/
size_t hash_key(VolumeDescriptor* key);

}//	end of namespace

#endif
