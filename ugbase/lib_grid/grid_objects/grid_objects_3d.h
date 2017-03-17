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

#ifndef __H__UG__grid_objects_3d__
#define __H__UG__grid_objects_3d__

#include "../grid/grid.h"
#include "common/math/ugmath.h"
#include "common/assert.h"
#include "grid_objects_0d.h"
#include "grid_objects_1d.h"
#include "grid_objects_2d.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
///	These numbers define where in the volume-section-container a volume will be stored.
/**	The order of the constants must not be changed! Algorithms may exist that rely on it.*/
enum VolumeContainerSections
{
	CSVOL_NONE = -1,
	CSVOL_TETRAHEDRON = 0,
	CSVOL_HEXAHEDRON = 1,
	CSVOL_PRISM = 2,
	CSVOL_PYRAMID = 3,
	CSVOL_OCTAHEDRON = 4
};

////////////////////////////////////////////////////////////////////////
//	TetrahedronDescriptor
///	only used to initialize a tetrahedron. for all other tasks you should use VolumeDescripor.
/**
 * please be sure to pass the vertices in the correct order:
 * v1, v2, v3: bottom-vertices in counterclockwise order (if viewed from the top).
 * v4: top
 */
class UG_API TetrahedronDescriptor
{
	public:
		TetrahedronDescriptor()	{}
		TetrahedronDescriptor(const TetrahedronDescriptor& td);
		TetrahedronDescriptor(const VolumeVertices& vv);
		TetrahedronDescriptor(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4);

		inline uint num_vertices() const	{return 4;}
		inline Vertex* vertex(size_t index) const	{return m_vertex[index];}

	protected:
		Vertex*	m_vertex[4];
};

////////////////////////////////////////////////////////////////////////
//	Tetrahedron
///	the most simple volume-element.
/**
 * order of vertices should be the same as described in \sa TetrahedronDescriptor
 *
 * \ingroup lib_grid_grid_objects
 */
class UG_API Tetrahedron : public Volume
{
	public:
		typedef Volume BaseClass;

		static const size_t NUM_VERTICES = 4;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<Tetrahedron*>(pObj) != NULL;}

		Tetrahedron()	{}
		Tetrahedron(const TetrahedronDescriptor& td);
		Tetrahedron(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4);

		virtual GridObject* create_empty_instance() const	{return new Tetrahedron;}

		virtual Vertex* vertex(size_t index) const	{return m_vertices[index];}
		virtual ConstVertexArray vertices() const		{return m_vertices;}
		virtual size_t num_vertices() const				{return 4;}

		virtual EdgeDescriptor edge_desc(int index) const;
		virtual void edge_desc(int index, EdgeDescriptor& edOut) const;
		virtual uint num_edges() const;

		virtual FaceDescriptor face_desc(int index) const;
		virtual void face_desc(int index, FaceDescriptor& fdOut) const;
		virtual uint num_faces() const;

		virtual Edge* create_edge(int index);	///< create the edge with index i and return it.
		virtual Face* create_face(int index);		///< create the face with index i and return it.

		virtual void get_vertex_indices_of_edge(size_t& ind1Out,
												size_t& ind2Out,
												size_t edgeInd) const;

		virtual void get_vertex_indices_of_face(std::vector<size_t>& indsOut,
												size_t side) const;

		virtual int get_edge_index_from_vertices(	const size_t vi0,
													const size_t vi1) const;

		virtual int get_face_edge_index (	const size_t faceInd,
											const size_t faceEdgeInd) const;

		virtual std::pair<GridBaseObjectId, int> get_opposing_object(Vertex* vrt) const;

	///	Creates new volume elements through refinement.
	/**	Make sure that newEdgeVertices contains 6 vertex pointers.
	 *	newFaceVertices is ignored for Tetrahedrons.*/
		virtual bool refine(std::vector<Volume*>& vNewVolumesOut,
							Vertex** ppNewVertexOut,
							Vertex** newEdgeVertices,
							Vertex** newFaceVertices,
							Vertex* newVolumeVertex,
							const Vertex& prototypeVertex,
							Vertex** pSubstituteVertices = NULL,
							vector3* corners = NULL,
							bool* isSnapPoint = NULL);

		virtual bool is_regular_ref_rule(int edgeMarks) const;

		virtual bool collapse_edge(std::vector<Volume*>& vNewVolumesOut,
								int edgeIndex, Vertex* newVertex,
								std::vector<Vertex*>* pvSubstituteVertices = NULL);

		virtual void get_flipped_orientation(VolumeDescriptor& vdOut) const;

		virtual int container_section() const
		{return CSVOL_TETRAHEDRON;}

		virtual ReferenceObjectID reference_object_id() const
		{return ROID_TETRAHEDRON;}

	protected:
		virtual void set_vertex(uint index, Vertex* pVrt)	{m_vertices[index] = pVrt;}

	protected:
		Vertex*	m_vertices[4];
};

template <>
class geometry_traits<Tetrahedron>
{
	public:
		typedef GenericGridObjectIterator<Tetrahedron*, VolumeIterator>		iterator;
		typedef ConstGenericGridObjectIterator<Tetrahedron*, VolumeIterator,
															ConstVolumeIterator>	const_iterator;

		typedef TetrahedronDescriptor Descriptor;
		typedef Volume 		grid_base_object;

		enum
		{
			CONTAINER_SECTION = CSVOL_TETRAHEDRON,
			BASE_OBJECT_ID = VOLUME
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_TETRAHEDRON;
};

typedef geometry_traits<Tetrahedron>::iterator			TetrahedronIterator;
typedef geometry_traits<Tetrahedron>::const_iterator	ConstTetrahedronIterator;



////////////////////////////////////////////////////////////////////////
//	HexahedronDescriptor
///	only used to initialize a hexahedron. for all other tasks you should use VolumeDescripor.
/**
 * please be sure to pass the vertices in the correct order:
 * v1, v2, v3, v4: bottom-vertices in counterclockwise order (if viewed from the top).
 * v5, v6, v7, v8: top-vertices in counterclockwise order (if viewed from the top).
 */
class UG_API HexahedronDescriptor
{
	public:
		HexahedronDescriptor()	{}
		HexahedronDescriptor(const HexahedronDescriptor& td);
		HexahedronDescriptor(const VolumeVertices& vv);
		HexahedronDescriptor(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4,
							Vertex* v5, Vertex* v6, Vertex* v7, Vertex* v8);

		inline uint num_vertices() const	{return 8;}
		inline Vertex* vertex(size_t index) const	{return m_vertex[index];}

	protected:
		Vertex*	m_vertex[8];
};

////////////////////////////////////////////////////////////////////////
//	Hexahedron
///	A volume element with 6 quadrilateral sides.
/**
 * Order of vertices should be the same as described in \sa HexahedronDescriptor
 *
 * \ingroup lib_grid_grid_objects
 */
class UG_API Hexahedron : public Volume
{
	public:
		typedef Volume BaseClass;

		static const size_t NUM_VERTICES = 8;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<Hexahedron*>(pObj) != NULL;}

		Hexahedron()	{}
		Hexahedron(const HexahedronDescriptor& td);
		Hexahedron(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4,
					Vertex* v5, Vertex* v6, Vertex* v7, Vertex* v8);

		virtual GridObject* create_empty_instance() const	{return new Hexahedron;}

		virtual Vertex* vertex(size_t index) const	{return m_vertices[index];}
		virtual ConstVertexArray vertices() const		{return m_vertices;}
		virtual size_t num_vertices() const				{return 8;}

		virtual EdgeDescriptor edge_desc(int index) const;
		virtual void edge_desc(int index, EdgeDescriptor& edOut) const;
		virtual uint num_edges() const;

		virtual FaceDescriptor face_desc(int index) const;
		virtual void face_desc(int index, FaceDescriptor& fdOut) const;
		virtual uint num_faces() const;

		virtual Edge* create_edge(int index);	///< create the edge with index i and return it.
		virtual Face* create_face(int index);		///< create the face with index i and return it.

		virtual void get_vertex_indices_of_edge(size_t& ind1Out,
												size_t& ind2Out,
												size_t edgeInd) const;
		
		virtual void get_vertex_indices_of_face(std::vector<size_t>& indsOut,
												size_t side) const;

		virtual int get_edge_index_from_vertices(	const size_t vi0,
													const size_t vi1) const;

		virtual int get_face_edge_index (	const size_t faceInd,
											const size_t faceEdgeInd) const;

		virtual bool get_opposing_side(FaceVertices* f, FaceDescriptor& fdOut) const;

		virtual std::pair<GridBaseObjectId, int> get_opposing_object(Vertex* vrt) const;

	///	see Volume::refine for a detailed description.
		virtual bool refine(std::vector<Volume*>& vNewVolumesOut,
							Vertex** ppNewVertexOut,
							Vertex** newEdgeVertices,
							Vertex** newFaceVertices,
							Vertex* newVolumeVertex,
							const Vertex& prototypeVertex,
							Vertex** pSubstituteVertices = NULL,
							vector3* corners = NULL,
							bool* isSnapPoint = NULL);

		virtual bool is_regular_ref_rule(int edgeMarks) const;

		virtual bool collapse_edge(std::vector<Volume*>& vNewVolumesOut,
								int edgeIndex, Vertex* newVertex,
								std::vector<Vertex*>* pvSubstituteVertices = NULL);

		virtual void get_flipped_orientation(VolumeDescriptor& vdOut) const;

		virtual int container_section() const	{return CSVOL_HEXAHEDRON;}
		virtual ReferenceObjectID reference_object_id() const {return ROID_HEXAHEDRON;}

	protected:
		virtual void set_vertex(uint index, Vertex* pVrt)	{m_vertices[index] = pVrt;}

	protected:
		Vertex*	m_vertices[8];
};

template <>
class geometry_traits<Hexahedron>
{
	public:
		typedef GenericGridObjectIterator<Hexahedron*, VolumeIterator>			iterator;
		typedef ConstGenericGridObjectIterator<Hexahedron*, VolumeIterator,
															 ConstVolumeIterator>	const_iterator;

		typedef HexahedronDescriptor Descriptor;
		typedef Volume 		grid_base_object;

		enum
		{
			CONTAINER_SECTION = CSVOL_HEXAHEDRON,
			BASE_OBJECT_ID = VOLUME
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_HEXAHEDRON;
};

typedef geometry_traits<Hexahedron>::iterator	HexahedronIterator;
typedef geometry_traits<Hexahedron>::const_iterator	ConstHexahedronIterator;


////////////////////////////////////////////////////////////////////////
//	PrismDescriptor
///	only used to initialize a prism. for all other tasks you should use VolumeDescripor.
/**
 * please be sure to pass the vertices in the correct order:
 * v1, v2, v3: bottom-vertices in counterclockwise order (if viewed from the top).
 * v4, v5, v6: top-vertices in counterclockwise order (if viewed from the top).
 */
class UG_API PrismDescriptor
{
	public:
		PrismDescriptor()	{}
		PrismDescriptor(const PrismDescriptor& td);
		PrismDescriptor(const VolumeVertices& vv);
		PrismDescriptor(Vertex* v1, Vertex* v2, Vertex* v3,
						Vertex* v4, Vertex* v5, Vertex* v6);

		inline uint num_vertices() const	{return 6;}
		inline Vertex* vertex(size_t index) const	{return m_vertex[index];}

	protected:
		Vertex*	m_vertex[6];
};

////////////////////////////////////////////////////////////////////////
//	Prism
///	A volume element with 2 triangle and 3 quadrilateral sides.
/**
 * order of vertices should be the same as described in \sa PrismDescriptor
 *
 * \ingroup lib_grid_grid_objects
 */
class UG_API Prism : public Volume
{
	public:
		typedef Volume BaseClass;

		static const size_t NUM_VERTICES = 6;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<Prism*>(pObj) != NULL;}

		Prism()	{}
		Prism(const PrismDescriptor& td);
		Prism(Vertex* v1, Vertex* v2, Vertex* v3,
				Vertex* v4, Vertex* v5, Vertex* v6);

		virtual GridObject* create_empty_instance() const	{return new Prism;}

		virtual Vertex* vertex(size_t index) const	{return m_vertices[index];}
		virtual ConstVertexArray vertices() const		{return m_vertices;}
		virtual size_t num_vertices() const				{return 6;}

		virtual EdgeDescriptor edge_desc(int index) const;
		virtual void edge_desc(int index, EdgeDescriptor& edOut) const;
		virtual uint num_edges() const;

		virtual FaceDescriptor face_desc(int index) const;
		virtual void face_desc(int index, FaceDescriptor& fdOut) const;
		virtual uint num_faces() const;

		virtual Edge* create_edge(int index);	///< create the edge with index i and return it.
		virtual Face* create_face(int index);		///< create the face with index i and return it.

		virtual void get_vertex_indices_of_edge(size_t& ind1Out,
												size_t& ind2Out,
												size_t edgeInd) const;

		virtual void get_vertex_indices_of_face(std::vector<size_t>& indsOut,
												size_t side) const;

		virtual int get_edge_index_from_vertices(	const size_t vi0,
													const size_t vi1) const;

		virtual int get_face_edge_index (	const size_t faceInd,
											const size_t faceEdgeInd) const;

		virtual bool get_opposing_side(FaceVertices* f, FaceDescriptor& fdOut) const;

		virtual std::pair<GridBaseObjectId, int> get_opposing_object(Vertex* vrt) const;

	///	see Volume::refine for a detailed description.
		virtual bool refine(std::vector<Volume*>& vNewVolumesOut,
							Vertex** ppNewVertexOut,
							Vertex** newEdgeVertices,
							Vertex** newFaceVertices,
							Vertex* newVolumeVertex,
							const Vertex& prototypeVertex,
							Vertex** pSubstituteVertices = NULL,
							vector3* corners = NULL,
							bool* isSnapPoint = NULL);

		virtual bool is_regular_ref_rule(int edgeMarks) const;

		virtual bool collapse_edge(std::vector<Volume*>& vNewVolumesOut,
								int edgeIndex, Vertex* newVertex,
								std::vector<Vertex*>* pvSubstituteVertices = NULL);

		virtual void get_flipped_orientation(VolumeDescriptor& vdOut) const;

		virtual int container_section() const	{return CSVOL_PRISM;}
		virtual ReferenceObjectID reference_object_id() const {return ROID_PRISM;}

	protected:
		virtual void set_vertex(uint index, Vertex* pVrt)	{m_vertices[index] = pVrt;}

	protected:
		Vertex*	m_vertices[6];
};

template <>
class geometry_traits<Prism>
{
	public:
		typedef GenericGridObjectIterator<Prism*, VolumeIterator>				iterator;
		typedef ConstGenericGridObjectIterator<Prism*, VolumeIterator,
															 ConstVolumeIterator>	const_iterator;

		typedef PrismDescriptor Descriptor;
		typedef Volume 		grid_base_object;

		enum
		{
			CONTAINER_SECTION = CSVOL_PRISM,
			BASE_OBJECT_ID = VOLUME
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_PRISM;
};

typedef geometry_traits<Prism>::iterator		PrismIterator;
typedef geometry_traits<Prism>::const_iterator	ConstPrismIterator;


////////////////////////////////////////////////////////////////////////
//	PyramidDescriptor
///	only used to initialize a pyramids. for all other tasks you should use VolumeDescripor.
/**
 * please be sure to pass the vertices in the correct order:
 * v1, v2, v3, v4: bottom-vertices in counterclockwise order (if viewed from the top).
 * v5: top-vertex.
 */
class UG_API PyramidDescriptor
{
	public:
		PyramidDescriptor()	{}
		PyramidDescriptor(const PyramidDescriptor& td);
		PyramidDescriptor(const VolumeVertices& vv);
		PyramidDescriptor(Vertex* v1, Vertex* v2, Vertex* v3,
						Vertex* v4, Vertex* v5);

		inline uint num_vertices() const	{return 5;}
		inline Vertex* vertex(size_t index) const	{return m_vertex[index];}

	protected:
		Vertex*	m_vertex[5];
};

////////////////////////////////////////////////////////////////////////
//	Pyramid
///	A volume element with 4 triangle and 1 quadrilateral sides.
/**
 * order of vertices should be the same as described in \sa PyramidDescriptor
 *
 * \ingroup lib_grid_grid_objects
 */
class UG_API Pyramid : public Volume
{
	public:
		typedef Volume BaseClass;

		static const size_t NUM_VERTICES = 5;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<Pyramid*>(pObj) != NULL;}

		Pyramid()	{}
		Pyramid(const PyramidDescriptor& td);
		Pyramid(Vertex* v1, Vertex* v2, Vertex* v3,
				Vertex* v4, Vertex* v5);

		virtual GridObject* create_empty_instance() const	{return new Pyramid;}

		virtual Vertex* vertex(size_t index) const	{return m_vertices[index];}
		virtual ConstVertexArray vertices() const		{return m_vertices;}
		virtual size_t num_vertices() const				{return 5;}

		virtual EdgeDescriptor edge_desc(int index) const;
		virtual void edge_desc(int index, EdgeDescriptor& edOut) const;
		virtual uint num_edges() const;

		virtual FaceDescriptor face_desc(int index) const;
		virtual void face_desc(int index, FaceDescriptor& fdOut) const;
		virtual uint num_faces() const;

		virtual Edge* create_edge(int index);	///< create the edge with index i and return it.
		virtual Face* create_face(int index);		///< create the face with index i and return it.

		virtual void get_vertex_indices_of_edge(size_t& ind1Out,
												size_t& ind2Out,
												size_t edgeInd) const;

		virtual void get_vertex_indices_of_face(std::vector<size_t>& indsOut,
												size_t side) const;

		virtual int get_edge_index_from_vertices(	const size_t vi0,
													const size_t vi1) const;

		virtual int get_face_edge_index (	const size_t faceInd,
											const size_t faceEdgeInd) const;

		virtual std::pair<GridBaseObjectId, int> get_opposing_object(Vertex* vrt) const;

	///	see Volume::refine for a detailed description.
		virtual bool refine(std::vector<Volume*>& vNewVolumesOut,
							Vertex** ppNewVertexOut,
							Vertex** newEdgeVertices,
							Vertex** newFaceVertices,
							Vertex* newVolumeVertex,
							const Vertex& prototypeVertex,
							Vertex** pSubstituteVertices = NULL,
							vector3* corners = NULL,
							bool* isSnapPoint = NULL);

		virtual bool is_regular_ref_rule(int edgeMarks) const;

		virtual bool collapse_edge(std::vector<Volume*>& vNewVolumesOut,
								int edgeIndex, Vertex* newVertex,
								std::vector<Vertex*>* pvSubstituteVertices = NULL);

		virtual void get_flipped_orientation(VolumeDescriptor& vdOut) const;

		virtual int container_section() const	{return CSVOL_PYRAMID;}
		virtual ReferenceObjectID reference_object_id() const {return ROID_PYRAMID;}

	protected:
		virtual void set_vertex(uint index, Vertex* pVrt)	{m_vertices[index] = pVrt;}

	protected:
		Vertex*	m_vertices[5];
};

template <>
class geometry_traits<Pyramid>
{
	public:
		typedef GenericGridObjectIterator<Pyramid*, VolumeIterator>			iterator;
		typedef ConstGenericGridObjectIterator<Pyramid*, VolumeIterator,
															 ConstVolumeIterator>	const_iterator;

		typedef PyramidDescriptor Descriptor;
		typedef Volume 		grid_base_object;

		enum
		{
			CONTAINER_SECTION = CSVOL_PYRAMID,
			BASE_OBJECT_ID = VOLUME
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_PYRAMID;
};

typedef geometry_traits<Pyramid>::iterator			PyramidIterator;
typedef geometry_traits<Pyramid>::const_iterator	ConstPyramidIterator;



////////////////////////////////////////////////////////////////////////
//	OctahedronDescriptor
///	only used to initialize a octahedron. for all other tasks you should use VolumeDescripor.
/**
 * please be sure to pass the vertices in the correct order:
 * v1: bottom-vertex
 * v2, v3, v4, v5: middle-section-vertices in counterclockwise order (if viewed from the top).
 * v6: top-vertex
 */
class UG_API OctahedronDescriptor
{
	public:
		OctahedronDescriptor()	{}
		OctahedronDescriptor(const OctahedronDescriptor& td);
		OctahedronDescriptor(const VolumeVertices& vv);
		OctahedronDescriptor(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4, Vertex* v5, Vertex* v6);

		inline uint num_vertices() const	{return 6;}
		inline Vertex* vertex(size_t index) const	{return m_vertex[index];}

	protected:
		Vertex*	m_vertex[6];
};


////////////////////////////////////////////////////////////////////////
//	Octahedron
///	platonic solid with eight faces.
/**
 * order of vertices should be the same as described in \sa OctahedronDescriptor
 *
 * \ingroup lib_grid_grid_objects
 */
class UG_API Octahedron : public Volume
{
	public:
		typedef Volume BaseClass;

		static const size_t NUM_VERTICES = 6;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<Octahedron*>(pObj) != NULL;}

		Octahedron()	{}
		Octahedron(const OctahedronDescriptor& td);
		Octahedron(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4, Vertex* v5, Vertex* v6);

		virtual GridObject* create_empty_instance() const	{return new Octahedron;}

		virtual Vertex* vertex(size_t index) const	{return m_vertices[index];}
		virtual ConstVertexArray vertices() const		{return m_vertices;}
		virtual size_t num_vertices() const				{return 6;}

		virtual EdgeDescriptor edge_desc(int index) const;
		virtual void edge_desc(int index, EdgeDescriptor& edOut) const;
		virtual uint num_edges() const;

		virtual FaceDescriptor face_desc(int index) const;
		virtual void face_desc(int index, FaceDescriptor& fdOut) const;
		virtual uint num_faces() const;

		virtual Edge* create_edge(int index);	///< create the edge with index i and return it.
		virtual Face* create_face(int index);		///< create the face with index i and return it.

		virtual void get_vertex_indices_of_edge(size_t& ind1Out,
												size_t& ind2Out,
												size_t edgeInd) const;

		virtual void get_vertex_indices_of_face(std::vector<size_t>& indsOut,
												size_t side) const;

		virtual int get_edge_index_from_vertices(	const size_t vi0,
													const size_t vi1) const;

		virtual int get_face_edge_index (	const size_t faceInd,
											const size_t faceEdgeInd) const;

		virtual std::pair<GridBaseObjectId, int> get_opposing_object(Vertex* vrt) const;

	///	Creates new volume elements through refinement.
	/**	Make sure that newEdgeVertices contains 6 vertex pointers.
	 *	newFaceVertices is ignored for Octahedrons.*/
		virtual bool refine(std::vector<Volume*>& vNewVolumesOut,
							Vertex** ppNewVertexOut,
							Vertex** newEdgeVertices,
							Vertex** newFaceVertices,
							Vertex* newVolumeVertex,
							const Vertex& prototypeVertex,
							Vertex** pSubstituteVertices = NULL,
							vector3* corners = NULL,
							bool* isSnapPoint = NULL);

		virtual bool is_regular_ref_rule(int edgeMarks) const;

		virtual bool collapse_edge(std::vector<Volume*>& vNewVolumesOut,
								int edgeIndex, Vertex* newVertex,
								std::vector<Vertex*>* pvSubstituteVertices = NULL);

		virtual void get_flipped_orientation(VolumeDescriptor& vdOut) const;

		virtual int container_section() const	{return CSVOL_OCTAHEDRON;}
		virtual ReferenceObjectID reference_object_id() const {return ROID_OCTAHEDRON;}

	protected:
		virtual void set_vertex(uint index, Vertex* pVrt)	{m_vertices[index] = pVrt;}

	protected:
		Vertex*	m_vertices[6];
};

template <>
class geometry_traits<Octahedron>
{
	public:
		typedef GenericGridObjectIterator<Octahedron*, VolumeIterator>		iterator;
		typedef ConstGenericGridObjectIterator<Octahedron*, VolumeIterator,
															ConstVolumeIterator>	const_iterator;

		typedef OctahedronDescriptor Descriptor;
		typedef Volume 		grid_base_object;

		enum
		{
			CONTAINER_SECTION = CSVOL_OCTAHEDRON,
			BASE_OBJECT_ID = VOLUME
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_OCTAHEDRON;
};

typedef geometry_traits<Octahedron>::iterator			OctahedronIterator;
typedef geometry_traits<Octahedron>::const_iterator	ConstOctahedronIterator;

}//	end of namespace

#endif
