//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m09 d10

#ifndef __H__LIB_GRID__MULTI_GRID__
#define __H__LIB_GRID__MULTI_GRID__

#include <vector>
#include <cassert>
#include "grid/grid.h"
#include "tools/subset_handler_grid.h"
#include "common_attachments.h"
#include "common/util/array_util.h"

namespace ug
{
 
//	Predeclaration
class MultiGrid;

//	The following constants define the maximum number of children
//	for each element-type.
//	This makes sense, since it allows us to align all the element
//	information in memory.
const int MG_EDGE_MAX_EDGE_CHILDREN = 2;///< maximal number of edges that can be children of an edge.
const int MG_FACE_MAX_EDGE_CHILDREN = 4;///< maximal number of edges that can be children of a face.
const int MG_FACE_MAX_FACE_CHILDREN = 4;///< maximal number of faces that can be children of a face.
//	the values below are just guesses. They may be too high or to low.
//	They should be replaced by a more flexible system!
//const int MG_VOLUME_MAX_EDGE_CHILDREN = 32;///< maximal number of edges that can be children of a volume.
//const int MG_VOLUME_MAX_FACE_CHILDREN = 32;///< maximal number of faces that can be children of a volume.
//const int MG_VOLUME_MAX_VOLUME_CHILDREN = 10;///< maximal number of volumes that can be children of a volume.


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
	MGVertexInfo()		{clear();}
	inline void clear()	{m_pParent = m_pVrtChild = NULL;}
	inline bool has_children() const {return m_pVrtChild != NULL;}
	inline void add_child(VertexBase* elem)	{assert(!m_pVrtChild); m_pVrtChild = elem;}
	inline void remove_child(VertexBase* elem)	{m_pVrtChild = NULL;}
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
	MGEdgeInfo()		{clear();}
	inline void clear()	{m_pParent = m_pVrtChild = NULL; m_numEdgeChildren = 0;}
	inline bool has_children() const {return m_pVrtChild || m_numEdgeChildren;}
	inline void add_child(VertexBase* elem)	{assert(!m_pVrtChild); m_pVrtChild = elem;}
	inline void add_child(EdgeBase* elem)	{assert(m_numEdgeChildren < MG_EDGE_MAX_EDGE_CHILDREN); m_pEdgeChild[m_numEdgeChildren++] = elem;}
	inline void remove_child(VertexBase* elem)	{m_pVrtChild = NULL;}
	inline void remove_child(EdgeBase* elem)	{m_numEdgeChildren = ArrayEraseEntry(m_pEdgeChild, elem, m_numEdgeChildren);}
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
	byte				m_numEdgeChildren;///< primarily required during refinement
};

///	Holds information about face relations. Used internally.
/**	No parent included, since MGFaceInfos are not stored for surface elements.*/
struct MGFaceInfo
{
	MGFaceInfo()		{clear();}
	inline void clear()	{m_pVrtChild = NULL; m_numEdgeChildren = m_numFaceChildren = 0;}
	inline bool has_children() const {return m_pVrtChild || m_numEdgeChildren || m_numFaceChildren;}
	inline void add_child(VertexBase* elem)	{assert(!m_pVrtChild); m_pVrtChild = elem;}
	inline void add_child(EdgeBase* elem)	{assert(m_numEdgeChildren < MG_FACE_MAX_EDGE_CHILDREN); m_pEdgeChild[m_numEdgeChildren++] = elem;}
	inline void add_child(Face* elem)		{assert(m_numFaceChildren < MG_FACE_MAX_FACE_CHILDREN); m_pFaceChild[m_numFaceChildren++] = elem;}
	inline void remove_child(VertexBase* elem)	{m_pVrtChild = NULL;}
	inline void remove_child(EdgeBase* elem)	{m_numEdgeChildren = ArrayEraseEntry(m_pEdgeChild, elem, m_numEdgeChildren);}
	inline void remove_child(Face* elem)		{m_numFaceChildren = ArrayEraseEntry(m_pFaceChild, elem, m_numFaceChildren);}
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
	byte				m_numEdgeChildren;///< primarily required during refinement
	byte				m_numFaceChildren;///< primarily required during refinement
};

///	Holds information about volume relations. Used internally.
/**	No parent included, since MGFaceInfos are not stored for surface elements.*/
struct MGVolumeInfo
{
	MGVolumeInfo() : m_pVrtChild(NULL)		{}
	inline void clear()	{m_pVrtChild = NULL; m_edgeChildren.clear(); m_faceChildren.clear(); m_volumeChildren.clear();}
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
	//EdgeBase* 			m_pEdgeChild[MG_VOLUME_MAX_EDGE_CHILDREN];
	//Face*				m_pFaceChild[MG_VOLUME_MAX_FACE_CHILDREN];
	//Volume*				m_pVolChild[MG_VOLUME_MAX_VOLUME_CHILDREN];
	//byte				m_numEdgeChildren;///< primarily required during refinement
	//byte				m_numFaceChildren;///< primarily required during refinement
	//byte				m_numVolChildren;///< primarily required during refinement
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

////////////////////////////////////////////////////////////////////////
/**
 * Inherits from \sa Grid.
 *
 * A MultiGrid represents a grid hierarchy. Elements in a level have a
 * parent / child relationship to elements in lower / higher levels.
 * Such a hierarchy is normally created by repeated refinement of a
 * coarse grid.
 * Enhances the Grid interface by methods that work on specific levels.
 * The MultiGrid stores all elements in one grid.
 * The hierarchy is managed by a SubsetHandler.
 * If elements are created and hierarchical insertion is activated, then
 * new elements are added one layer higher than their parents.
 * (NULL indicates base-level).
 *
 * Whenever a a level is added or removed, a message is posted at the
 * associated MessageHub (MultiGrid::message_hub()). The message has the type
 * GridMessage_MultiGridChanged (defined in "lib_grid/lib_grid_messages.h").
 * You may register a callback at the grids message-hub if you want to react
 * on such a message.
 *
 * \ingroup lib_grid_grid
 */
class MultiGrid : public Grid, public GridObserver
{
	friend struct MGVertexInfo;
	friend struct MGEdgeInfo;
	friend struct MGFaceInfo;
	friend struct MGVolumeInfo;

	protected:
		typedef MGVertexInfo VertexInfo;
		typedef MGEdgeInfo EdgeInfo;
		typedef MGFaceInfo FaceInfo;
		typedef MGVolumeInfo VolumeInfo;

	public:
		//	methods from Grid, that would be hidden if not explicitly
		//	declared as required.
		using Grid::begin;
		using Grid::end;
		using Grid::num;
		using Grid::get_geometric_objects;
		using Grid::create;
		using Grid::create_by_cloning;

	public:
		MultiGrid();
	///	initialises the grid with the given option.
	/**	pass an or-combination of constants enumerated in
	 *  VertexOptions, EdgeOptions, FaceOptions, VolumeOptions and GridOptions.*/
		MultiGrid(uint options);
		
		virtual ~MultiGrid();

		void enable_hierarchical_insertion(bool bEnable);
		inline bool hierarchical_insertion_enabled()	{return m_bHierarchicalInsertion;}

	////////////////////////////////////////////////
	//	element creation
	///	create a custom element on a specific level.
	/**
	 * TGeomObj has to be a geometric object type as described in geometric_base_objects.h.
	 * This method should only be used if a geometric object has to be created
	 * without a parent in higher levels of the hierarchy.
	 * Use the create method derived from ug::Grid if you want to specify a parent.
	 * \{ */
		template<class TGeomObj>
		typename geometry_traits<TGeomObj>::iterator
		create(size_t level);

		template <class TGeomObj>
		typename geometry_traits<TGeomObj>::iterator
		create(const typename geometry_traits<TGeomObj>::Descriptor& descriptor,
				size_t level);
	/**	\} */

	///	this method creates a new vertex, which has the same type as pCloneMe.
		VertexBaseIterator create_by_cloning(VertexBase* pCloneMe,
											int level);

	///	this method creates a new edge, which has the same type as pCloneMe.
		EdgeBaseIterator create_by_cloning(EdgeBase* pCloneMe,
										   const EdgeVertices& ev,
										   int level);

	///	this method creates a new face, which has the same type as pCloneMe.
		FaceIterator create_by_cloning(Face* pCloneMe,
									   const FaceVertices& fv,
									   int level);

	///	this method creates a new volume, which has the same type as pCloneMe.
		VolumeIterator create_by_cloning(Volume* pCloneMe,
										 const VolumeVertices& vv,
										 int level);

	///	number of levels
		inline size_t num_levels() const	{return (size_t)m_hierarchy.num_subsets();}

	///	index of the highest level.
		inline size_t top_level() const;

	///	creates new (empty) levels until num_levels() == lvl+1
		inline void level_required(int lvl);

		template <class TElem> inline
		size_t num(int level) const		{return m_hierarchy.num<TElem>(level);}

		template <class TElem> inline
		typename geometry_traits<TElem>::iterator
		begin(int level)
		{
			//assert(level < (int)num_levels() && "ERROR in MultiGrid::begin(...): requested level too high!");
			if(level >= (int)num_levels())
				return end<TElem>();
			return m_hierarchy.begin<TElem>(level);
		}

		template <class TElem> inline
		typename geometry_traits<TElem>::iterator
		end(int level)
		{
			//assert(level < (int)num_levels() && "ERROR in MultiGrid::end(...): requested level too high!");
			if(level >= (int)num_levels())
				return end<TElem>();
			return m_hierarchy.end<TElem>(level);
		}

		template <class TElem> inline
		typename geometry_traits<TElem>::const_iterator
		begin(int level) const
		{
			//assert(level < (int)num_levels() && "ERROR in MultiGrid::begin(...): requested level too high!");
			if(level >= (int)num_levels())
				return end<TElem>();
			return m_hierarchy.begin<TElem>(level);
		}

		template <class TElem> inline
		typename geometry_traits<TElem>::const_iterator
		end(int level) const
		{
			//assert(level < (int)num_levels() && "ERROR in MultiGrid::end(...): requested level too high!");
			if(level >= (int)num_levels())
				return end<TElem>();
			return m_hierarchy.end<TElem>(level);
		}

	//	geometric-object-collection
		inline GeometricObjectCollection
		get_geometric_objects(int level)
		{return m_hierarchy.get_geometric_objects_in_subset(level);}
		
	//	multi-level-geometric-object-collection
		virtual GeometricObjectCollection get_geometric_objects()
		{return m_hierarchy.get_geometric_objects();}
		
		template <class TElem> inline
		int get_level(TElem* elem) const
		{return m_hierarchy.get_subset_index(elem);}

		GeometricObject* get_parent(GeometricObject* parent) const;
		inline GeometricObject* get_parent(VertexBase* o) const	{return get_info(o).m_pParent;}
		inline GeometricObject* get_parent(EdgeBase* o) const	{return get_info(o).m_pParent;}
		inline GeometricObject* get_parent(Face* o) const		{return m_aaParentFACE[o];}
		inline GeometricObject* get_parent(Volume* o) const		{return m_aaParentVOL[o];}

	//	number of children
		template <class TElem> inline
		bool has_children(TElem* elem) const
		{return get_info(elem).has_children();}

	////////////////////////////////
	//	CHILD QUANTITIES
	///	returns the number of children of the given child-type
	/** \{ */
		template <class TChild, class TElem>
		inline size_t num_children(TElem* elem)	const	{return num_children(elem, TChild());}

		template <class TChild>
		size_t num_children(GeometricObject* elem)	const;
	/** \} */

	///	returns the total number of children and grand-children.
	/**	Only children of the same type as the given elements are regarded here.*/
		template <class TElem>
		inline size_t num_children_total(TElem* elem)	const;

	///	Returns the number of child vertices
		template <class TElem>
		inline size_t num_child_vertices(TElem* elem) const	{return get_info(elem).num_child_vertices();}

	///	Returns the number of child edges
	/**	\{	*/
		template <class TElem>
		inline size_t num_child_edges(TElem* elem) const	{return get_info(elem).num_child_edges();}
		inline size_t num_child_edges(VertexBase*) const	{return 0;}
	/**	\}	*/

	///	Returns the number of child faces
	/**	\{	*/
		template <class TElem>
		inline size_t num_child_faces(TElem* elem) const	{return get_info(elem).num_child_faces();}
		inline size_t num_child_faces(VertexBase*) const	{return 0;}
		inline size_t num_child_faces(EdgeBase*) const		{return 0;}
	/**	\}	*/

	///	Returns the number of child volumes
	/**	\{	*/
		inline size_t num_child_volumes(Volume* elem) const	{return get_info(elem).num_child_volumes();}
		template <class TElem>
		inline size_t num_child_volumes(TElem*) const		{return 0;}
	/**	\}	*/


	////////////////////////////////
	//	CHILD ACCESS
	///	returns the i-th child of the given child-type
	/** \{ */
		template <class TChild, class TElem>
		inline TChild* get_child(TElem* elem, size_t ind) const	{return get_child(elem, ind, TChild());}

		template <class TChild>
		TChild* get_child(GeometricObject* elem, size_t ind) const;
	/** \} */

	///	Returns the child vertex of the given element or NULL if there is none
		template <class TElem>
		inline VertexBase* get_child_vertex(TElem* elem) const	{return get_info(elem).child_vertex();}

	///	Returns the child edges of the given element or NULL if there is none
	/**	\{	*/
		template <class TElem>
		inline EdgeBase* get_child_edge(TElem* elem, size_t ind) const	{return get_info(elem).child_edge(ind);}
		inline EdgeBase* get_child_edge(VertexBase*, size_t) const		{return NULL;}
	/**	\}	*/

	///	Returns the child faces of the given element or NULL if there is none
	/**	\{	*/
		template <class TElem>
		inline Face* get_child_face(TElem* elem, size_t ind) const	{return get_info(elem).child_face(ind);}
		inline Face* get_child_face(VertexBase*, size_t) const		{return NULL;}
		inline Face* get_child_face(EdgeBase*, size_t) const			{return NULL;}
	/**	\}	*/

	///	Returns the child volumes of the given element or NULL if there is none
	/**	\{	*/
		inline Volume* get_child_volume(Volume* elem, size_t ind) const	{return get_info(elem).child_volume(ind);}
		template <class TElem>
		inline Volume* get_child_volume(TElem*, size_t) const	{return NULL;}
	/**	\}	*/

	///	for debug purposes
		void check_edge_elem_infos(int level) const;
	///	for debug purposes
		void check_face_elem_infos(int level) const;
	///	for debug purposes
		void check_volume_elem_infos(int level) const;
		
	///	this method may be removed in future versions of the MultiGrid-class.
	/**	You really shouldn't use this method!!!*/
		SubsetHandler& get_hierarchy_handler()		{return m_hierarchy;}
		
	////////////////////////////////////////////////////////////////////////
	//	Don't invoke the following methods directly!
	//	They are intended for internal feedback only.

	//	grid callbacks
		virtual void elements_to_be_cleared(Grid* grid);

	 /**  In order to correctly register vrt in the hierarchy, we have to
	 *  replace pParent with vrt in the list of children of pParents parent.
	 *  This means that if a grid-observer registered after the multi-grid itself,
	 *  iterates over the list of children of pParents parent in its
	 *  vertex_created method, it won't find pParent.*/
		virtual void vertex_created(Grid* grid, VertexBase* vrt,
									GeometricObject* pParent = NULL,
									bool replacesParent = false);

	 /**  In order to correctly register e in the hierarchy, we have to
	 *  replace pParent with e in the list of children of pParents parent.
	 *  This means that if a grid-observer registered after the multi-grid itself,
	 *  iterates over the list of children of pParents parent in its
	 *  edge_created method, it won't find pParent.*/
		virtual void edge_created(Grid* grid, EdgeBase* e,
									GeometricObject* pParent = NULL,
									bool replacesParent = false);

	 /**  In order to correctly register f in the hierarchy, we have to
	 *  replace pParent with f in the list of children of pParents parent.
	 *  This means that if a grid-observer registered after the multi-grid itself,
	 *  iterates over the list of children of pParents parent in its
	 *  face_created method, it won't find pParent.*/
		virtual void face_created(Grid* grid, Face* f,
									GeometricObject* pParent = NULL,
									bool replacesParent = false);

	 /**  In order to correctly register vol in the hierarchy, we have to
	 *  replace pParent with vol in the list of children of pParents parent.
	 *  This means that if a grid-observer registered after the multi-grid itself,
	 *  iterates over the list of children of pParents parent in its
	 *  volume_created method, it won't find pParent.*/
		virtual void volume_created(Grid* grid, Volume* vol,
									GeometricObject* pParent = NULL,
									bool replacesParent = false);

		virtual void vertex_to_be_erased(Grid* grid, VertexBase* vrt,
										 VertexBase* replacedBy = NULL);

		virtual void edge_to_be_erased(Grid* grid, EdgeBase* e,
										 EdgeBase* replacedBy = NULL);

		virtual void face_to_be_erased(Grid* grid, Face* f,
										 Face* replacedBy = NULL);

		virtual void volume_to_be_erased(Grid* grid, Volume* vol,
										 Volume* replacedBy = NULL);

	protected:

	//	Note: VertexInfo and EdgeInfo are stored directly, FaceInfo and
	//	VolumeInfo are stored dynamically.
		typedef Attachment<GeometricObject*>	AParent;
		typedef Attachment<VertexInfo>			AVertexInfo;
		typedef Attachment<EdgeInfo>			AEdgeInfo;
		typedef Attachment<FaceInfo*>			AFaceInfo;
		typedef Attachment<VolumeInfo*>			AVolumeInfo;

	protected:
	//	initialization
		void init();
		
	//	create levels
		void create_levels(int numLevels);

	//	info-access
		inline VertexInfo& get_info(VertexBase* v);
		inline EdgeInfo& get_info(EdgeBase* e);
		inline FaceInfo& get_info(Face* f);
		inline VolumeInfo& get_info(Volume* v);

	//	const info-access
		inline const VertexInfo& get_info(VertexBase* v) const;
		inline const EdgeInfo& get_info(EdgeBase* e) const;
		inline const FaceInfo& get_info(Face* f) const;
		inline const VolumeInfo& get_info(Volume* v) const;

	//	elem creation
		template <class TElem>
		inline void element_created(TElem* elem)	{element_created<TElem, TElem>(elem, NULL);}

		template <class TElem, class TParent>
		void element_created(TElem* elem, TParent* pParent);

	///	called if a newly created element shall replace an old one
		template <class TElem, class TParent>
		void element_created(TElem* elem, TParent* pParent, TElem* pReplaceMe);

	///	this method is called for elements that havn't got any parent.
		template <class TElem>
		void element_to_be_erased(TElem* elem);

	///	this method is called for elements with a parent.
		template <class TElem, class TParent>
		void element_to_be_erased(TElem* elem, TParent* pParent);

		//template <class TElem>
		//void element_to_be_replaced(TElem* elemOld, TElem* elemNew);


	///	returning the number of children of the type of the dummy-argument.
	/**	\{ */
		template <class TElem>
		inline size_t num_children(TElem* elem, const VertexBase&) const
			{return num_child_vertices(elem);}

		template <class TElem>
		inline size_t num_children(TElem* elem, const EdgeBase&) const
			{return num_child_edges(elem);}

		template <class TElem>
		inline size_t num_children(TElem* elem, const Face&) const
			{return num_child_faces(elem);}

		template <class TElem>
		inline size_t num_children(TElem* elem, const Volume&) const
			{return num_child_volumes(elem);}
	/**	\} */

	///	returning the i-th child of the type of the dummy-argument.
	/**	\{ */
		template <class TElem>
		inline VertexBase* get_child(TElem* elem, size_t ind, const VertexBase&) const
			{return get_child_vertex(elem);}

		template <class TElem>
		inline EdgeBase* get_child(TElem* elem, size_t ind, const EdgeBase&) const
			{return get_child_edge(elem, ind);}

		template <class TElem>
		inline Face* get_child(TElem* elem, size_t ind, const Face&) const
			{return get_child_face(elem, ind);}

		template <class TElem>
		inline Volume* get_child(TElem* elem, size_t ind, const Volume&) const
			{return get_child_volume(elem, ind);}
	/**	\} */

	///	sets the parent for the given object
	/**	\{ */
		inline void set_parent(VertexBase* o, GeometricObject* p)	{get_info(o).m_pParent = p;}
		inline void set_parent(EdgeBase* o, GeometricObject* p)		{get_info(o).m_pParent = p;}
		inline void set_parent(Face* o, GeometricObject* p)			{m_aaParentFACE[o] = p;}
		inline void set_parent(Volume* o, GeometricObject* p)		{m_aaParentVOL[o] = p;}
	/**	\} */

	///	creates the info-object for the given object (if necessary)
	/**	\{ */
		inline void create_child_info(VertexBase* o){}
		inline void create_child_info(EdgeBase* o)	{}
		inline void create_child_info(Face* o)		{if(!m_aaFaceInf[o]) m_aaFaceInf[o] = new FaceInfo();}
		inline void create_child_info(Volume* o)	{if(!m_aaVolInf[o]) m_aaVolInf[o] = new VolumeInfo();}
	/**	\} */

	///	releases the info-object for the given object (if necessary)
	/**	\{ */
		inline void release_child_info(VertexBase* o)	{}
		inline void release_child_info(EdgeBase* o)		{}
		inline void release_child_info(Face* o)			{if(m_aaFaceInf[o]) delete m_aaFaceInf[o]; m_aaFaceInf[o] = NULL;}
		inline void release_child_info(Volume* o)		{if(m_aaVolInf[o]) delete m_aaVolInf[o]; m_aaVolInf[o] = NULL;}
	/**	\} */

	protected:
	//	hierarchy
		SubsetHandler	m_hierarchy;
		bool m_bHierarchicalInsertion;

	//	message-id for changed-messages
		int			m_msgId;

	//	parent attachment
		AParent		m_aParent;

	//	info attachments
		AVertexInfo	m_aVertexInfo;
		AEdgeInfo	m_aEdgeInfo;
		AFaceInfo	m_aFaceInfo;
		AVolumeInfo	m_aVolumeInfo;

	//	parent access - only required for faces and volumes.
		Grid::FaceAttachmentAccessor<AParent>		m_aaParentFACE;
		Grid::VolumeAttachmentAccessor<AParent>		m_aaParentVOL;

	//	element info access
		Grid::VertexAttachmentAccessor<AVertexInfo>	m_aaVrtInf;
		Grid::EdgeAttachmentAccessor<AEdgeInfo>		m_aaEdgeInf;
		Grid::FaceAttachmentAccessor<AFaceInfo>		m_aaFaceInf;
		Grid::VolumeAttachmentAccessor<AVolumeInfo>	m_aaVolInf;
};



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///	wrapper that allows to write method that can operate on a ug::Grid and a ug::MultiGrid.
/**
 * This template class is specialized for ug::Grid and ug::MultiGrid.
 *
 * The MGWrapper will most likely be enhanced with more mg-methods, like
 * has_child, get_child or parent.
 */
template <class TGrid>
class MGWrapper
{
	public:
		MGWrapper(TGrid& grid);
		
		inline uint num_levels() const;

		template <class TElem> inline
		uint num(int level) const;

		template <class TElem> inline
		typename geometry_traits<TElem>::iterator
		begin(int level);

		template <class TElem> inline
		typename geometry_traits<TElem>::iterator
		end(int level);
};

}//	end of namespace

////////////////////////////////
//	include implementation
#include "multi_grid_impl.hpp"

#endif
