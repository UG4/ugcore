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

#include "grid_objects_0d.h"
#include "grid_objects_1d.h"
#include "grid_objects_2d.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
///	These numbers define where in the volume-section-container a volume will be stored.
/**	The order of the constants must not be changed! Algorithms may exist that rely on it.*/
enum VolumeContainerSectionsw
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
		TetrahedronDescriptor()	= default;
		TetrahedronDescriptor(const TetrahedronDescriptor& td);
		explicit TetrahedronDescriptor(const VolumeVertices& vv);
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
		using BaseClass = Volume;

		static constexpr size_t NUM_VERTICES = 4;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<Tetrahedron*>(pObj) != nullptr;}

		Tetrahedron() = default;
		explicit Tetrahedron(const TetrahedronDescriptor& td);
		Tetrahedron(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4);

		[[nodiscard]] GridObject* create_empty_instance() const override {return new Tetrahedron;}

		[[nodiscard]] Vertex* vertex(size_t index) const override {return m_vertices[index];}
		[[nodiscard]] ConstVertexArray vertices() const override {return m_vertices;}
		[[nodiscard]] size_t num_vertices() const override {return 4;}

		[[nodiscard]] EdgeDescriptor edge_desc(int index) const override;

		void edge_desc(int index, EdgeDescriptor& edOut) const override;

		[[nodiscard]] uint num_edges() const override;

		[[nodiscard]] FaceDescriptor face_desc(int index) const override;

		void face_desc(int index, FaceDescriptor& fdOut) const override;

		[[nodiscard]] uint num_faces() const override;

		Edge* create_edge(int index) override;	///< create the edge with index i and return it.
		Face* create_face(int index) override;		///< create the face with index i and return it.

		void get_vertex_indices_of_edge(size_t& ind1Out,
		                                size_t& ind2Out,
		                                size_t edgeInd) const override;

		void get_vertex_indices_of_face(std::vector<size_t>& indsOut,
		                                size_t side) const override;

		[[nodiscard]] int get_edge_index_from_vertices(	const size_t vi0,
			                                 const size_t vi1) const override;

		[[nodiscard]] int get_face_edge_index (	const size_t faceInd,
			                         const size_t faceEdgeInd) const override;

		std::pair<GridBaseObjectId, int> get_opposing_object(Vertex* vrt) const override;

	///	Creates new volume elements through refinement.
	/**	Make sure that newEdgeVertices contains 6 vertex pointers.
	 *	newFaceVertices is ignored for Tetrahedrons.*/
		bool refine(std::vector<Volume*>& vNewVolumesOut,
		            Vertex** ppNewVertexOut,
		            Vertex** newEdgeVertices,
		            Vertex** newFaceVertices,
		            Vertex* newVolumeVertex,
		            const Vertex& prototypeVertex,
		            Vertex** pSubstituteVertices = nullptr,
		            vector3* corners = nullptr,
		            bool* isSnapPoint = nullptr) override;

		[[nodiscard]] bool is_regular_ref_rule(int edgeMarks) const override;

		bool collapse_edge(std::vector<Volume*>& vNewVolumesOut,
		                   int edgeIndex, Vertex* newVertex,
		                   std::vector<Vertex*>* pvSubstituteVertices = nullptr) override;

		void get_flipped_orientation(VolumeDescriptor& vdOut) const override;

		[[nodiscard]] int container_section() const override {return VolumeContainerSectionsw::CSVOL_TETRAHEDRON;}

		[[nodiscard]] ReferenceObjectID reference_object_id() const override {return ReferenceObjectID::ROID_TETRAHEDRON;}

	protected:
		void set_vertex(uint index, Vertex* pVrt) override {m_vertices[index] = pVrt;}

	protected:
		Vertex*	m_vertices[4];
};

template <>
class geometry_traits<Tetrahedron>
{
	public:
		using iterator = GenericGridObjectIterator<Tetrahedron*, VolumeIterator>;
		using const_iterator = ConstGenericGridObjectIterator<Tetrahedron*, VolumeIterator,
			ConstVolumeIterator>;

		using Descriptor = TetrahedronDescriptor;
		using grid_base_object = Volume;

		enum
		{
			CONTAINER_SECTION = VolumeContainerSectionsw::CSVOL_TETRAHEDRON,
			BASE_OBJECT_ID = GridBaseObjectId::VOLUME
		};
		static constexpr ReferenceObjectID REFERENCE_OBJECT_ID = ReferenceObjectID::ROID_TETRAHEDRON;
};

using TetrahedronIterator = geometry_traits<Tetrahedron>::iterator;
using ConstTetrahedronIterator = geometry_traits<Tetrahedron>::const_iterator;



////////////////////////////////////////////////////////////////////////
//	HexahedronDescriptor
///	only used to initialize a hexahedron. for all other tasks you should use VolumeDescriptor.
/**
 * please be sure to pass the vertices in the correct order:
 * v1, v2, v3, v4: bottom-vertices in counterclockwise order (if viewed from the top).
 * v5, v6, v7, v8: top-vertices in counterclockwise order (if viewed from the top).
 */
class UG_API HexahedronDescriptor
{
	public:
		HexahedronDescriptor() = default;
		HexahedronDescriptor(const HexahedronDescriptor& td);
		explicit HexahedronDescriptor(const VolumeVertices& vv);
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
		using BaseClass = Volume;

		static constexpr size_t NUM_VERTICES = 8;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<Hexahedron*>(pObj) != nullptr;}

		Hexahedron() = default;
		explicit Hexahedron(const HexahedronDescriptor& td);
		Hexahedron(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4,
					Vertex* v5, Vertex* v6, Vertex* v7, Vertex* v8);

		[[nodiscard]] GridObject* create_empty_instance() const override {return new Hexahedron;}

		[[nodiscard]] Vertex* vertex(size_t index) const override {return m_vertices[index];}
		[[nodiscard]] ConstVertexArray vertices() const override {return m_vertices;}
		[[nodiscard]] size_t num_vertices() const override {return 8;}

		[[nodiscard]] EdgeDescriptor edge_desc(int index) const override;

		void edge_desc(int index, EdgeDescriptor& edOut) const override;

		[[nodiscard]] uint num_edges() const override;

		[[nodiscard]] FaceDescriptor face_desc(int index) const override;

		void face_desc(int index, FaceDescriptor& fdOut) const override;

		[[nodiscard]] uint num_faces() const override;

		Edge* create_edge(int index) override;	///< create the edge with index i and return it.
		Face* create_face(int index) override;		///< create the face with index i and return it.

		void get_vertex_indices_of_edge(size_t& ind1Out,
		                                size_t& ind2Out,
		                                size_t edgeInd) const override;

		void get_vertex_indices_of_face(std::vector<size_t>& indsOut,
		                                size_t side) const override;

		[[nodiscard]] int get_edge_index_from_vertices(	const size_t vi0,
			                                 const size_t vi1) const override;

		[[nodiscard]] int get_face_edge_index (	const size_t faceInd,
			                         const size_t faceEdgeInd) const override;

		bool get_opposing_side(FaceVertices* f, FaceDescriptor& fdOut) const override;

		std::pair<GridBaseObjectId, int> get_opposing_object(Vertex* vrt) const override;

	///	see Volume::refine for a detailed description.
		bool refine(std::vector<Volume*>& vNewVolumesOut,
		            Vertex** ppNewVertexOut,
		            Vertex** newEdgeVertices,
		            Vertex** newFaceVertices,
		            Vertex* newVolumeVertex,
		            const Vertex& prototypeVertex,
		            Vertex** pSubstituteVertices = nullptr,
		            vector3* corners = nullptr,
		            bool* isSnapPoint = nullptr) override;

		[[nodiscard]] bool is_regular_ref_rule(int edgeMarks) const override;

		bool collapse_edge(std::vector<Volume*>& vNewVolumesOut,
		                   int edgeIndex, Vertex* newVertex,
		                   std::vector<Vertex*>* pvSubstituteVertices = nullptr) override;

		void get_flipped_orientation(VolumeDescriptor& vdOut) const override;

		[[nodiscard]] int container_section() const override {return VolumeContainerSectionsw::CSVOL_HEXAHEDRON;}
		[[nodiscard]] ReferenceObjectID reference_object_id() const override {return ReferenceObjectID::ROID_HEXAHEDRON;}

	protected:
		void set_vertex(uint index, Vertex* pVrt) override {m_vertices[index] = pVrt;}

	protected:
		Vertex*	m_vertices[8];
};

template <>
class geometry_traits<Hexahedron>
{
	public:
		using iterator = GenericGridObjectIterator<Hexahedron*, VolumeIterator>;
		using const_iterator = ConstGenericGridObjectIterator<Hexahedron*, VolumeIterator,
			ConstVolumeIterator>;

		using Descriptor = HexahedronDescriptor;
		using grid_base_object = Volume;

		enum
		{
			CONTAINER_SECTION = VolumeContainerSectionsw::CSVOL_HEXAHEDRON,
			BASE_OBJECT_ID = GridBaseObjectId::VOLUME
		};
		static constexpr ReferenceObjectID REFERENCE_OBJECT_ID = ReferenceObjectID::ROID_HEXAHEDRON;
};

using HexahedronIterator = geometry_traits<Hexahedron>::iterator;
using ConstHexahedronIterator = geometry_traits<Hexahedron>::const_iterator;


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
		PrismDescriptor() = default;
		PrismDescriptor(const PrismDescriptor& td);
		explicit PrismDescriptor(const VolumeVertices& vv);
		PrismDescriptor(Vertex* v1, Vertex* v2, Vertex* v3,
						Vertex* v4, Vertex* v5, Vertex* v6);

		inline uint num_vertices() const {return 6;}
		inline Vertex* vertex(size_t index) const {return m_vertex[index];}

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
		using BaseClass = Volume;

		static constexpr size_t NUM_VERTICES = 6;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<Prism*>(pObj) != nullptr;}

		Prism()	= default;
		explicit Prism(const PrismDescriptor& td);
		Prism(Vertex* v1, Vertex* v2, Vertex* v3,
				Vertex* v4, Vertex* v5, Vertex* v6);

		[[nodiscard]] GridObject* create_empty_instance() const override {return new Prism;}

		[[nodiscard]] Vertex* vertex(size_t index) const override {return m_vertices[index];}
		[[nodiscard]] ConstVertexArray vertices() const override {return m_vertices;}
		[[nodiscard]] size_t num_vertices() const override {return 6;}

		[[nodiscard]] EdgeDescriptor edge_desc(int index) const override;

		void edge_desc(int index, EdgeDescriptor& edOut) const override;

		[[nodiscard]] uint num_edges() const override;

		[[nodiscard]] FaceDescriptor face_desc(int index) const override;

		void face_desc(int index, FaceDescriptor& fdOut) const override;

		[[nodiscard]] uint num_faces() const override;

		Edge* create_edge(int index) override;	///< create the edge with index i and return it.
		Face* create_face(int index) override;		///< create the face with index i and return it.

		void get_vertex_indices_of_edge(size_t& ind1Out,
		                                size_t& ind2Out,
		                                size_t edgeInd) const override;

		void get_vertex_indices_of_face(std::vector<size_t>& indsOut,
		                                size_t side) const override;

		[[nodiscard]] int get_edge_index_from_vertices(	const size_t vi0,
			                                 const size_t vi1) const override;

		[[nodiscard]] int get_face_edge_index (	const size_t faceInd,
			                         const size_t faceEdgeInd) const override;

		bool get_opposing_side(FaceVertices* f, FaceDescriptor& fdOut) const override;

		std::pair<GridBaseObjectId, int> get_opposing_object(Vertex* vrt) const override;

	///	see Volume::refine for a detailed description.
		bool refine(std::vector<Volume*>& vNewVolumesOut,
		            Vertex** ppNewVertexOut,
		            Vertex** newEdgeVertices,
		            Vertex** newFaceVertices,
		            Vertex* newVolumeVertex,
		            const Vertex& prototypeVertex,
		            Vertex** pSubstituteVertices = nullptr,
		            vector3* corners = nullptr,
		            bool* isSnapPoint = nullptr) override;

		[[nodiscard]] bool is_regular_ref_rule(int edgeMarks) const override;

		bool collapse_edge(std::vector<Volume*>& vNewVolumesOut,
		                   int edgeIndex, Vertex* newVertex,
		                   std::vector<Vertex*>* pvSubstituteVertices = nullptr) override;

		void get_flipped_orientation(VolumeDescriptor& vdOut) const override;

		[[nodiscard]] int container_section() const override {return VolumeContainerSectionsw::CSVOL_PRISM;}
		[[nodiscard]] ReferenceObjectID reference_object_id() const override {return ReferenceObjectID::ROID_PRISM;}

	protected:
		void set_vertex(uint index, Vertex* pVrt) override {m_vertices[index] = pVrt;}

	protected:
		Vertex*	m_vertices[6];
};

template <>
class geometry_traits<Prism>
{
	public:
		using iterator = GenericGridObjectIterator<Prism*, VolumeIterator>;
		using const_iterator = ConstGenericGridObjectIterator<Prism*, VolumeIterator,
			ConstVolumeIterator>;

		using Descriptor = PrismDescriptor;
		using grid_base_object = Volume;

		enum
		{
			CONTAINER_SECTION = VolumeContainerSectionsw::CSVOL_PRISM,
			BASE_OBJECT_ID = GridBaseObjectId::VOLUME
		};
		static constexpr ReferenceObjectID REFERENCE_OBJECT_ID = ReferenceObjectID::ROID_PRISM;
};

using PrismIterator = geometry_traits<Prism>::iterator;
using ConstPrismIterator = geometry_traits<Prism>::const_iterator;


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
		PyramidDescriptor()	= default;
		PyramidDescriptor(const PyramidDescriptor& td);
		explicit PyramidDescriptor(const VolumeVertices& vv);
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
		using BaseClass = Volume;

		static constexpr size_t NUM_VERTICES = 5;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<Pyramid*>(pObj) != nullptr;}

		Pyramid() = default;
		explicit Pyramid(const PyramidDescriptor& td);
		Pyramid(Vertex* v1, Vertex* v2, Vertex* v3,
				Vertex* v4, Vertex* v5);

		[[nodiscard]] GridObject* create_empty_instance() const override {return new Pyramid;}

		[[nodiscard]] Vertex* vertex(size_t index) const override {return m_vertices[index];}
		[[nodiscard]] ConstVertexArray vertices() const override {return m_vertices;}
		[[nodiscard]] size_t num_vertices() const override {return 5;}

		[[nodiscard]] EdgeDescriptor edge_desc(int index) const override;

		void edge_desc(int index, EdgeDescriptor& edOut) const override;

		[[nodiscard]] uint num_edges() const override;

		[[nodiscard]] FaceDescriptor face_desc(int index) const override;

		void face_desc(int index, FaceDescriptor& fdOut) const override;

		[[nodiscard]] uint num_faces() const override;

		Edge* create_edge(int index) override;	///< create the edge with index i and return it.
		Face* create_face(int index) override;		///< create the face with index i and return it.

		void get_vertex_indices_of_edge(size_t& ind1Out,
		                                size_t& ind2Out,
		                                size_t edgeInd) const override;

		void get_vertex_indices_of_face(std::vector<size_t>& indsOut,
		                                size_t side) const override;

		[[nodiscard]] int get_edge_index_from_vertices(	const size_t vi0,
			                                 const size_t vi1) const override;

		[[nodiscard]] int get_face_edge_index (	const size_t faceInd,
			                         const size_t faceEdgeInd) const override;

		std::pair<GridBaseObjectId, int> get_opposing_object(Vertex* vrt) const override;

	///	see Volume::refine for a detailed description.
		bool refine(std::vector<Volume*>& vNewVolumesOut,
		            Vertex** ppNewVertexOut,
		            Vertex** newEdgeVertices,
		            Vertex** newFaceVertices,
		            Vertex* newVolumeVertex,
		            const Vertex& prototypeVertex,
		            Vertex** pSubstituteVertices = nullptr,
		            vector3* corners = nullptr,
		            bool* isSnapPoint = nullptr) override;

		[[nodiscard]] bool is_regular_ref_rule(int edgeMarks) const override;

		bool collapse_edge(std::vector<Volume*>& vNewVolumesOut,
		                   int edgeIndex, Vertex* newVertex,
		                   std::vector<Vertex*>* pvSubstituteVertices = nullptr) override;

		void get_flipped_orientation(VolumeDescriptor& vdOut) const override;

		[[nodiscard]] int container_section() const override {return VolumeContainerSectionsw::CSVOL_PYRAMID;}
		[[nodiscard]] ReferenceObjectID reference_object_id() const override {return ReferenceObjectID::ROID_PYRAMID;}

	protected:
		void set_vertex(uint index, Vertex* pVrt) override {m_vertices[index] = pVrt;}

	protected:
		Vertex*	m_vertices[5];
};

template <>
class geometry_traits<Pyramid>
{
	public:
		using iterator = GenericGridObjectIterator<Pyramid*, VolumeIterator>;
		using const_iterator = ConstGenericGridObjectIterator<Pyramid*, VolumeIterator,
			ConstVolumeIterator>;

		using Descriptor = PyramidDescriptor;
		using grid_base_object = Volume;

		enum
		{
			CONTAINER_SECTION = VolumeContainerSectionsw::CSVOL_PYRAMID,
			BASE_OBJECT_ID = GridBaseObjectId::VOLUME
		};
		static constexpr ReferenceObjectID REFERENCE_OBJECT_ID = ReferenceObjectID::ROID_PYRAMID;
};

using PyramidIterator = geometry_traits<Pyramid>::iterator;
using ConstPyramidIterator = geometry_traits<Pyramid>::const_iterator;



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
		OctahedronDescriptor() = default;
		OctahedronDescriptor(const OctahedronDescriptor& td);
		explicit OctahedronDescriptor(const VolumeVertices& vv);
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
		using BaseClass = Volume;

		static constexpr size_t NUM_VERTICES = 6;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<Octahedron*>(pObj) != nullptr;}

		Octahedron() = default;

		explicit Octahedron(const OctahedronDescriptor& td);
		Octahedron(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4, Vertex* v5, Vertex* v6);

		[[nodiscard]] GridObject* create_empty_instance() const override {return new Octahedron;}

		[[nodiscard]] Vertex* vertex(size_t index) const override {return m_vertices[index];}
		[[nodiscard]] ConstVertexArray vertices() const override {return m_vertices;}
		[[nodiscard]] size_t num_vertices() const override {return 6;}

		[[nodiscard]] EdgeDescriptor edge_desc(int index) const override;

		void edge_desc(int index, EdgeDescriptor& edOut) const override;

		[[nodiscard]] uint num_edges() const override;

		[[nodiscard]] FaceDescriptor face_desc(int index) const override;

		void face_desc(int index, FaceDescriptor& fdOut) const override;

		[[nodiscard]] uint num_faces() const override;

		Edge* create_edge(int index) override;	///< create the edge with index i and return it.
		Face* create_face(int index) override;		///< create the face with index i and return it.

		void get_vertex_indices_of_edge(size_t& ind1Out,
		                                size_t& ind2Out,
		                                size_t edgeInd) const override;

		void get_vertex_indices_of_face(std::vector<size_t>& indsOut,
		                                size_t side) const override;

		[[nodiscard]] int get_edge_index_from_vertices(	const size_t vi0,
			                                 const size_t vi1) const override;

		[[nodiscard]] int get_face_edge_index (	const size_t faceInd,
			                         const size_t faceEdgeInd) const override;

		std::pair<GridBaseObjectId, int> get_opposing_object(Vertex* vrt) const override;

	///	Creates new volume elements through refinement.
	/**	Make sure that newEdgeVertices contains 6 vertex pointers.
	 *	newFaceVertices is ignored for Octahedrons.*/
		bool refine(std::vector<Volume*>& vNewVolumesOut,
		            Vertex** ppNewVertexOut,
		            Vertex** newEdgeVertices,
		            Vertex** newFaceVertices,
		            Vertex* newVolumeVertex,
		            const Vertex& prototypeVertex,
		            Vertex** pSubstituteVertices = nullptr,
		            vector3* corners = nullptr,
		            bool* isSnapPoint = nullptr) override;

		[[nodiscard]] bool is_regular_ref_rule(int edgeMarks) const override;

		bool collapse_edge(std::vector<Volume*>& vNewVolumesOut,
		                   int edgeIndex, Vertex* newVertex,
		                   std::vector<Vertex*>* pvSubstituteVertices = nullptr) override;

		void get_flipped_orientation(VolumeDescriptor& vdOut) const override;

		[[nodiscard]] int container_section() const override {return VolumeContainerSectionsw::CSVOL_OCTAHEDRON;}
		[[nodiscard]] ReferenceObjectID reference_object_id() const override {return ReferenceObjectID::ROID_OCTAHEDRON;}

	protected:
		void set_vertex(uint index, Vertex* pVrt) override {m_vertices[index] = pVrt;}

	protected:
		Vertex*	m_vertices[6];
};

template <>
class geometry_traits<Octahedron>
{
	public:
		using iterator = GenericGridObjectIterator<Octahedron*, VolumeIterator>;
		using const_iterator = ConstGenericGridObjectIterator<Octahedron*, VolumeIterator,
			ConstVolumeIterator>;

		using Descriptor = OctahedronDescriptor;
		using grid_base_object = Volume;

		enum
		{
			CONTAINER_SECTION = VolumeContainerSectionsw::CSVOL_OCTAHEDRON,
			BASE_OBJECT_ID = GridBaseObjectId::VOLUME
		};
		static constexpr ReferenceObjectID REFERENCE_OBJECT_ID = ReferenceObjectID::ROID_OCTAHEDRON;
};

using OctahedronIterator = geometry_traits<Octahedron>::iterator;
using ConstOctahedronIterator = geometry_traits<Octahedron>::const_iterator;

}//	end of namespace

#endif
