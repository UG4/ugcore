/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * UML_LOOK
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

#ifndef __H__UG__multi_grid_child_info__
#define __H__UG__multi_grid_child_info__

#include <vector>
#include "common/error.h"
#include "common/assert.h"
#include "lib_grid/grid/grid_base_objects.h"
namespace ug{


//	Predeclaration
class MultiGrid;

//	The following constants define the maximum number of children
//	for edges amd faces.
//	Note that no constants are defined for volumes since the range of
//	possible children is too high (especially during adaptive refinement)
constexpr int MG_EDGE_MAX_EDGE_CHILDREN = 2;///< maximal number of edges that can be children of an edge.
constexpr int MG_FACE_MAX_EDGE_CHILDREN = 4;///< maximal number of edges that can be children of a face.
constexpr int MG_FACE_MAX_FACE_CHILDREN = 4;///< maximal number of faces that can be children of a face.

///	Holds information about vertex relations. Used internally.
/**
 * a vertex can have one of the following states:
 * MGES_NORMAL, MGES_CONSTRAINED, MGES_FIXED
 * the following should always hold true:
 * children of a constrained node are constrained, too.
 * children of a fixed node are fixed, too.
 */
struct MGVertexInfo
{
	MGVertexInfo() : m_pParent(nullptr), m_pVrtChild(nullptr)		{}
	inline void clear_children() {m_pVrtChild = nullptr;}
	inline bool has_children() const {return m_pVrtChild != nullptr;}
	inline void add_child(Vertex* elem)	 {assert(!m_pVrtChild); m_pVrtChild = elem;}
	inline void add_child(Edge*) {UG_THROW("INVALID OPERATION!");}//dummy...
	inline void add_child(Face*) {UG_THROW("INVALID OPERATION!");}//dummy...
	inline void add_child(Volume*) {UG_THROW("INVALID OPERATION!");}//dummy...
	inline void remove_child(Vertex* elem) {m_pVrtChild = nullptr;}
	inline void remove_child(Edge*) {UG_THROW("INVALID OPERATION!");}//dummy...
	inline void remove_child(Face*) {UG_THROW("INVALID OPERATION!");}//dummy...
	inline void remove_child(Volume*) {UG_THROW("INVALID OPERATION!");}//dummy...
	inline void replace_child(Vertex* elem, Vertex* child)	{assert(child == m_pVrtChild); m_pVrtChild = elem;}
	void unregister_from_children(MultiGrid& mg);
	[[nodiscard]] size_t num_child_vertices() const	{return m_pVrtChild ? 1 : 0;}

	[[nodiscard]] Vertex* child_vertex() const {return m_pVrtChild;}

	GridObject*	m_pParent;
private:
	Vertex*			m_pVrtChild;
};

///	Holds information about edge relations. Used internally.
/**
 * an edge can have one of the following states:
 * MGES_NORMAL, MGES_CONSTRAINING, MGES_CONSTRAINED
 * the following should always hold true:
 * children of a constraining edge are constrained.
 * children of a constrained edge are constrained.
 */
struct MGEdgeInfo
{
	MGEdgeInfo() : m_pParent(nullptr), m_pVrtChild(nullptr), m_numEdgeChildren(0)		{}
	inline void clear_children() {m_pVrtChild = nullptr; m_numEdgeChildren = 0;}
	[[nodiscard]] inline bool has_children() const {return m_pVrtChild || m_numEdgeChildren;}
	inline void add_child(Vertex* elem) {assert(!m_pVrtChild); m_pVrtChild = elem;}
	inline void add_child(Edge* elem) {assert(m_numEdgeChildren < MG_EDGE_MAX_EDGE_CHILDREN); m_pEdgeChild[m_numEdgeChildren++] = elem;}
	inline void add_child(Face*) {UG_THROW("INVALID OPERATION!");}//dummy...
	inline void add_child(Volume*) {UG_THROW("INVALID OPERATION!");}//dummy...
	inline void remove_child(Vertex* elem) {m_pVrtChild = nullptr;}
	inline void remove_child(Edge* elem) {m_numEdgeChildren = ArrayEraseEntry(m_pEdgeChild, elem, m_numEdgeChildren);}
	inline void remove_child(Face*) {UG_THROW("INVALID OPERATION!");}//dummy...
	inline void remove_child(Volume*) {UG_THROW("INVALID OPERATION!");}//dummy...
	inline void replace_child(Vertex* elem, Vertex* child) {assert(child == m_pVrtChild); m_pVrtChild = elem;}
	inline void replace_child(Edge* elem, Edge* child) {ArrayReplaceEntry(m_pEdgeChild, elem, child, m_numEdgeChildren);}
	void unregister_from_children(MultiGrid& mg);
	[[nodiscard]] size_t num_child_vertices() const {return m_pVrtChild ? 1 : 0;}
	[[nodiscard]] size_t num_child_edges() const {return m_numEdgeChildren;}

	[[nodiscard]] Vertex* child_vertex() const {return m_pVrtChild;}
	[[nodiscard]] Edge* child_edge(size_t i) const {assert(i < num_child_edges()); return m_pEdgeChild[i];}

	GridObject*	m_pParent;
private:
	Vertex* m_pVrtChild;
	Edge* m_pEdgeChild[MG_EDGE_MAX_EDGE_CHILDREN];
	byte_t m_numEdgeChildren;
};

///	Holds information about face relations. Used internally.
/**	No parent included, since MGFaceInfos are not stored for surface elements.*/
struct MGFaceInfo
{
	MGFaceInfo() : m_pVrtChild(nullptr), m_numEdgeChildren(0), m_numFaceChildren(0)	{}
	inline void clear_children() {m_pVrtChild = nullptr; m_numEdgeChildren = m_numFaceChildren = 0;}
	[[nodiscard]] inline bool has_children() const {return m_pVrtChild || m_numEdgeChildren || m_numFaceChildren;}
	inline void add_child(Vertex* elem)	{assert(!m_pVrtChild); m_pVrtChild = elem;}
	inline void add_child(Edge* elem) {assert(m_numEdgeChildren < MG_FACE_MAX_EDGE_CHILDREN); m_pEdgeChild[m_numEdgeChildren++] = elem;}
	inline void add_child(Face* elem) {assert(m_numFaceChildren < MG_FACE_MAX_FACE_CHILDREN); m_pFaceChild[m_numFaceChildren++] = elem;}
	inline void add_child(Volume*) {UG_THROW("INVALID OPERATION!");}//dummy...
	inline void remove_child(Vertex* elem) {m_pVrtChild = nullptr;}
	inline void remove_child(Edge* elem) {m_numEdgeChildren = ArrayEraseEntry(m_pEdgeChild, elem, m_numEdgeChildren);}
	inline void remove_child(Face* elem) {m_numFaceChildren = ArrayEraseEntry(m_pFaceChild, elem, m_numFaceChildren);}
	inline void remove_child(Volume*) {UG_THROW("INVALID OPERATION!");}//dummy...
	inline void replace_child(Vertex* elem, Vertex* child)	{assert(child == m_pVrtChild); m_pVrtChild = elem;}
	inline void replace_child(Edge* elem, Edge* child) {ArrayReplaceEntry(m_pEdgeChild, elem, child, m_numEdgeChildren);}
	inline void replace_child(Face* elem, Face* child) {ArrayReplaceEntry(m_pFaceChild, elem, child, m_numFaceChildren);}
	void unregister_from_children(MultiGrid& mg);

	[[nodiscard]] size_t num_child_vertices() const {return m_pVrtChild ? 1 : 0;}
	[[nodiscard]] size_t num_child_edges() const {return m_numEdgeChildren;}
	[[nodiscard]] size_t num_child_faces() const {return m_numFaceChildren;}

	[[nodiscard]] Vertex* child_vertex() const {return m_pVrtChild;}
	[[nodiscard]] Edge* child_edge(size_t i) const {assert(i < num_child_edges()); return m_pEdgeChild[i];}
	[[nodiscard]] Face* child_face(size_t i) const {assert(i < num_child_faces()); return m_pFaceChild[i];}

private:
	Vertex* m_pVrtChild;
	Edge* m_pEdgeChild[MG_FACE_MAX_EDGE_CHILDREN];
	Face* m_pFaceChild[MG_FACE_MAX_FACE_CHILDREN];
	byte_t m_numEdgeChildren;
	byte_t m_numFaceChildren;
};

///	Holds information about volume relations. Used internally.
/**	No parent included, since MGFaceInfos are not stored for surface elements.*/
struct MGVolumeInfo
{
	MGVolumeInfo() : m_pVrtChild(nullptr) {}
	inline void clear_children() {m_pVrtChild = nullptr; m_edgeChildren.clear(); m_faceChildren.clear(); m_volumeChildren.clear();}
	[[nodiscard]] inline bool has_children() const {return m_pVrtChild || !(m_edgeChildren.empty() && m_faceChildren.empty() && m_volumeChildren.empty());}
	inline void add_child(Vertex* elem) {assert(!m_pVrtChild); m_pVrtChild = elem;}
	inline void add_child(Edge* elem) {m_edgeChildren.push_back(elem);}
	inline void add_child(Face* elem) {m_faceChildren.push_back(elem);}
	inline void add_child(Volume* elem) {m_volumeChildren.push_back(elem);}
	inline void remove_child(Vertex* elem) {m_pVrtChild = nullptr;}
	inline void remove_child(Edge* elem) {ArraySwapWithLast(&m_edgeChildren.front(), elem, m_edgeChildren.size()); m_edgeChildren.pop_back();}
	inline void remove_child(Face* elem) {ArraySwapWithLast(&m_faceChildren.front(), elem, m_faceChildren.size()); m_faceChildren.pop_back();}
	inline void remove_child(Volume* elem) {ArraySwapWithLast(&m_volumeChildren.front(), elem, m_volumeChildren.size()); m_volumeChildren.pop_back();}
	inline void replace_child(Vertex* elem, Vertex* child) {assert(child == m_pVrtChild); m_pVrtChild = elem;}
	inline void replace_child(Edge* elem, Edge* child) {ArrayReplaceEntry(&m_edgeChildren.front(), elem, child, m_edgeChildren.size());}
	inline void replace_child(Face* elem, Face* child) {ArrayReplaceEntry(&m_faceChildren.front(), elem, child, m_faceChildren.size());}
	inline void replace_child(Volume* elem, Volume* child) {ArrayReplaceEntry(&m_volumeChildren.front(), elem, child, m_volumeChildren.size());}
	void unregister_from_children(MultiGrid& mg);

	[[nodiscard]] size_t num_child_vertices() const {return m_pVrtChild ? 1 : 0;}
	[[nodiscard]] size_t num_child_edges() const {return m_edgeChildren.size();}
	[[nodiscard]] size_t num_child_faces() const {return m_faceChildren.size();}
	[[nodiscard]] size_t num_child_volumes() const {return m_volumeChildren.size();}

	[[nodiscard]] Vertex* child_vertex() const {return m_pVrtChild;}
	[[nodiscard]] Edge* child_edge(size_t i) const {assert(i < num_child_edges()); return m_edgeChildren[i];}
	[[nodiscard]] Face* child_face(size_t i) const {assert(i < num_child_faces()); return m_faceChildren[i];}
	[[nodiscard]] Volume* child_volume(size_t i) const {assert(i < num_child_volumes()); return m_volumeChildren[i];}

private:
	Vertex* m_pVrtChild;
	std::vector<Edge*> m_edgeChildren;
	std::vector<Face*> m_faceChildren;
	std::vector<Volume*> m_volumeChildren;
};

///	access to connected types. used internally
/**
 * has to contain:
 * 	- type definition info_type
 */
template <typename TElem> class mginfo_traits{};

///	vertex info traits. used internally.
template <> class mginfo_traits<Vertex>
{
	public:
		using info_type = MGVertexInfo;
};

///	edge info traits. used internally.
template <> class mginfo_traits<Edge>
{
	public:
		using info_type = MGEdgeInfo;
};

///	face info traits. used internally.
template <> class mginfo_traits<Face>
{
	public:
		using info_type = MGFaceInfo;
};

///	volume info traits. used internally.
template <> class mginfo_traits<Volume>
{
	public:
		using info_type = MGVolumeInfo;
};

}// end of namespace

#endif
