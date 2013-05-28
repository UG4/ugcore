// created by Sebastian Reiter
// s.b.reiter@gmail.com
// 11.01.2013 (d,m,y)

#ifndef __H__UG__multi_grid_child_info__
#define __H__UG__multi_grid_child_info__

#include <vector>
#include "common/error.h"
#include "common/assert.h"
#include "lib_grid/grid/geometric_base_objects.h"
namespace ug{


//	Predeclaration
class MultiGrid;

//	The following constants define the maximum number of children
//	for edges amd faces.
//	Note that no constants are defined for volumes since the range of
//	possible children is too high (especially during adaptive refinement)
const int MG_EDGE_MAX_EDGE_CHILDREN = 2;///< maximal number of edges that can be children of an edge.
const int MG_FACE_MAX_EDGE_CHILDREN = 4;///< maximal number of edges that can be children of a face.
const int MG_FACE_MAX_FACE_CHILDREN = 4;///< maximal number of faces that can be children of a face.

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
	MGVertexInfo() : m_pParent(NULL), m_pVrtChild(NULL)		{}
	inline void clear_children()	{m_pVrtChild = NULL;}
	inline bool has_children() const {return m_pVrtChild != NULL;}
	inline void add_child(VertexBase* elem)		{assert(!m_pVrtChild); m_pVrtChild = elem;}
	inline void add_child(EdgeBase*)			{UG_THROW("INVALID OPERATION!");}//dummy...
	inline void add_child(Face*)				{UG_THROW("INVALID OPERATION!");}//dummy...
	inline void add_child(Volume*)				{UG_THROW("INVALID OPERATION!");}//dummy...
	inline void remove_child(VertexBase* elem)	{m_pVrtChild = NULL;}
	inline void remove_child(EdgeBase*)			{UG_THROW("INVALID OPERATION!");}//dummy...
	inline void remove_child(Face*)				{UG_THROW("INVALID OPERATION!");}//dummy...
	inline void remove_child(Volume*)			{UG_THROW("INVALID OPERATION!");}//dummy...
	inline void replace_child(VertexBase* elem, VertexBase* child)	{assert(child == m_pVrtChild); m_pVrtChild = elem;}
	void unregister_from_children(MultiGrid& mg);
	size_t num_child_vertices() const	{return m_pVrtChild ? 1 : 0;}

	VertexBase* child_vertex() const {return m_pVrtChild;}

	GeometricObject*	m_pParent;
private:
	VertexBase*			m_pVrtChild;
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
	MGEdgeInfo() : m_pParent(NULL), m_pVrtChild(NULL), m_numEdgeChildren(0)		{}
	inline void clear_children()	{m_pVrtChild = NULL; m_numEdgeChildren = 0;}
	inline bool has_children() const {return m_pVrtChild || m_numEdgeChildren;}
	inline void add_child(VertexBase* elem)	{assert(!m_pVrtChild); m_pVrtChild = elem;}
	inline void add_child(EdgeBase* elem)	{assert(m_numEdgeChildren < MG_EDGE_MAX_EDGE_CHILDREN); m_pEdgeChild[m_numEdgeChildren++] = elem;}
	inline void add_child(Face*)			{UG_THROW("INVALID OPERATION!");}//dummy...
	inline void add_child(Volume*)			{UG_THROW("INVALID OPERATION!");}//dummy...
	inline void remove_child(VertexBase* elem)	{m_pVrtChild = NULL;}
	inline void remove_child(EdgeBase* elem)	{m_numEdgeChildren = ArrayEraseEntry(m_pEdgeChild, elem, m_numEdgeChildren);}
	inline void remove_child(Face*)				{UG_THROW("INVALID OPERATION!");}//dummy...
	inline void remove_child(Volume*)			{UG_THROW("INVALID OPERATION!");}//dummy...
	inline void replace_child(VertexBase* elem, VertexBase* child)	{assert(child == m_pVrtChild); m_pVrtChild = elem;}
	inline void replace_child(EdgeBase* elem, EdgeBase* child)		{ArrayReplaceEntry(m_pEdgeChild, elem, child, m_numEdgeChildren);}
	void unregister_from_children(MultiGrid& mg);
	size_t num_child_vertices() const	{return m_pVrtChild ? 1 : 0;}
	size_t num_child_edges() const		{return m_numEdgeChildren;}

	VertexBase* child_vertex() const		{return m_pVrtChild;}
	EdgeBase* child_edge(size_t i) const	{assert(i < num_child_edges()); return m_pEdgeChild[i];}

	GeometricObject*	m_pParent;
private:
	VertexBase*			m_pVrtChild;
	EdgeBase* 			m_pEdgeChild[MG_EDGE_MAX_EDGE_CHILDREN];
	byte				m_numEdgeChildren;
};

///	Holds information about face relations. Used internally.
/**	No parent included, since MGFaceInfos are not stored for surface elements.*/
struct MGFaceInfo
{
	MGFaceInfo() : m_pVrtChild(NULL), m_numEdgeChildren(0), m_numFaceChildren(0)	{}
	inline void clear_children()	{m_pVrtChild = NULL; m_numEdgeChildren = m_numFaceChildren = 0;}
	inline bool has_children() const {return m_pVrtChild || m_numEdgeChildren || m_numFaceChildren;}
	inline void add_child(VertexBase* elem)	{assert(!m_pVrtChild); m_pVrtChild = elem;}
	inline void add_child(EdgeBase* elem)	{assert(m_numEdgeChildren < MG_FACE_MAX_EDGE_CHILDREN); m_pEdgeChild[m_numEdgeChildren++] = elem;}
	inline void add_child(Face* elem)		{assert(m_numFaceChildren < MG_FACE_MAX_FACE_CHILDREN); m_pFaceChild[m_numFaceChildren++] = elem;}
	inline void add_child(Volume*)			{UG_THROW("INVALID OPERATION!");}//dummy...
	inline void remove_child(VertexBase* elem)	{m_pVrtChild = NULL;}
	inline void remove_child(EdgeBase* elem)	{m_numEdgeChildren = ArrayEraseEntry(m_pEdgeChild, elem, m_numEdgeChildren);}
	inline void remove_child(Face* elem)		{m_numFaceChildren = ArrayEraseEntry(m_pFaceChild, elem, m_numFaceChildren);}
	inline void remove_child(Volume*)			{UG_THROW("INVALID OPERATION!");}//dummy...
	inline void replace_child(VertexBase* elem, VertexBase* child)	{assert(child == m_pVrtChild); m_pVrtChild = elem;}
	inline void replace_child(EdgeBase* elem, EdgeBase* child)		{ArrayReplaceEntry(m_pEdgeChild, elem, child, m_numEdgeChildren);}
	inline void replace_child(Face* elem, Face* child)				{ArrayReplaceEntry(m_pFaceChild, elem, child, m_numFaceChildren);}
	void unregister_from_children(MultiGrid& mg);

	size_t num_child_vertices() const	{return m_pVrtChild ? 1 : 0;}
	size_t num_child_edges() const		{return m_numEdgeChildren;}
	size_t num_child_faces() const		{return m_numFaceChildren;}

	VertexBase* child_vertex() const {return m_pVrtChild;}
	EdgeBase* child_edge(size_t i) const	{assert(i < num_child_edges()); return m_pEdgeChild[i];}
	Face* child_face(size_t i) const		{assert(i < num_child_faces()); return m_pFaceChild[i];}

private:
	VertexBase*			m_pVrtChild;
	EdgeBase* 			m_pEdgeChild[MG_FACE_MAX_EDGE_CHILDREN];
	Face*				m_pFaceChild[MG_FACE_MAX_FACE_CHILDREN];
	byte				m_numEdgeChildren;
	byte				m_numFaceChildren;
};

///	Holds information about volume relations. Used internally.
/**	No parent included, since MGFaceInfos are not stored for surface elements.*/
struct MGVolumeInfo
{
	MGVolumeInfo() : m_pVrtChild(NULL)		{}
	inline void clear_children()	{m_pVrtChild = NULL; m_edgeChildren.clear(); m_faceChildren.clear(); m_volumeChildren.clear();}
	inline bool has_children()	const {return m_pVrtChild || !(m_edgeChildren.empty() && m_faceChildren.empty() && m_volumeChildren.empty());}
	inline void add_child(VertexBase* elem)	{assert(!m_pVrtChild); m_pVrtChild = elem;}
	inline void add_child(EdgeBase* elem)	{m_edgeChildren.push_back(elem);}
	inline void add_child(Face* elem)		{m_faceChildren.push_back(elem);}
	inline void add_child(Volume* elem)		{m_volumeChildren.push_back(elem);}
	inline void remove_child(VertexBase* elem)	{m_pVrtChild = NULL;}
	inline void remove_child(EdgeBase* elem)	{ArraySwapWithLast(&m_edgeChildren.front(), elem, m_edgeChildren.size()); m_edgeChildren.pop_back();}
	//{m_numEdgeChildren = ArrayEraseEntry(m_pEdgeChild, elem, m_numEdgeChildren);}
	inline void remove_child(Face* elem)		{ArraySwapWithLast(&m_faceChildren.front(), elem, m_faceChildren.size()); m_faceChildren.pop_back();}
	//{m_numFaceChildren = ArrayEraseEntry(m_pFaceChild, elem, m_numFaceChildren);}
	inline void remove_child(Volume* elem)		{ArraySwapWithLast(&m_volumeChildren.front(), elem, m_volumeChildren.size()); m_volumeChildren.pop_back();}
	//{m_numVolChildren = ArrayEraseEntry(m_pVolChild, elem, m_numVolChildren);}
	inline void replace_child(VertexBase* elem, VertexBase* child)	{assert(child == m_pVrtChild); m_pVrtChild = elem;}
	inline void replace_child(EdgeBase* elem, EdgeBase* child)		{ArrayReplaceEntry(&m_edgeChildren.front(), elem, child, m_edgeChildren.size());}
	inline void replace_child(Face* elem, Face* child)				{ArrayReplaceEntry(&m_faceChildren.front(), elem, child, m_faceChildren.size());}
	inline void replace_child(Volume* elem, Volume* child)			{ArrayReplaceEntry(&m_volumeChildren.front(), elem, child, m_volumeChildren.size());}
	void unregister_from_children(MultiGrid& mg);

	size_t num_child_vertices() const	{return m_pVrtChild ? 1 : 0;}
	size_t num_child_edges() const		{return m_edgeChildren.size();}
	size_t num_child_faces() const		{return m_faceChildren.size();}
	size_t num_child_volumes() const	{return m_volumeChildren.size();}

	VertexBase* child_vertex() const 		{return m_pVrtChild;}
	EdgeBase* child_edge(size_t i) const	{assert(i < num_child_edges()); return m_edgeChildren[i];}
	Face* child_face(size_t i) const		{assert(i < num_child_faces()); return m_faceChildren[i];}
	Volume* child_volume(size_t i) const	{assert(i < num_child_volumes()); return m_volumeChildren[i];}

private:
	VertexBase*				m_pVrtChild;
	std::vector<EdgeBase*>	m_edgeChildren;
	std::vector<Face*>		m_faceChildren;
	std::vector<Volume*>	m_volumeChildren;
};

///	access to connected types. used internally
/**
 * has to contain:
 * 	- typedef info_type
 */
template <class TElem> class mginfo_traits{};

///	vertex info traits. used internally.
template <> class mginfo_traits<VertexBase>
{
	public:
		typedef MGVertexInfo	info_type;
};

///	edge info traits. used internally.
template <> class mginfo_traits<EdgeBase>
{
	public:
		typedef MGEdgeInfo	info_type;
};

///	face info traits. used internally.
template <> class mginfo_traits<Face>
{
	public:
		typedef MGFaceInfo	info_type;
};

///	volume info traits. used internally.
template <> class mginfo_traits<Volume>
{
	public:
		typedef MGVolumeInfo	info_type;
};

}// end of namespace

#endif
