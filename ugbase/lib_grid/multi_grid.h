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

//TODO:	Improve implementation of child-handling.
//		Think about only allowing one child-vertex per element.
//		Do we really need to be able to add children of any type
//		to a vertex (think about the overhead of a section-container).

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
const int MG_VOLUME_MAX_EDGE_CHILDREN = 6;///< maximal number of edges that can be children of a volume.
const int MG_VOLUME_MAX_FACE_CHILDREN = 12;///< maximal number of faces that can be children of a volume.
const int MG_VOLUME_MAX_VOLUME_CHILDREN = 8;///< maximal number of volumes that can be children of a volume.

////////////////////////////////////////////////////////////////////////
///	constants that describe the state of an element
/**
 * If not all elements are marked for refinement, then there will be
 * a boundary of the marked area.
 * Vertices in the new level that have a parent vertex that lies on
 * this boundary are regarded as fixed vertices and have the state MGES_FIXED.
 * Vertices and edges in the new level that have a parent edge that lies
 * on the mark-boundary are regarded as constraied objects and have
 * the state MGES_CONSTRAINED.
 * If the parent edge of such an element is not a constrained edge itself,
 * it has the state MGES_CONSTRAINING.
 * All other elements have the state MGES_NORMAL.
 */
enum MGElementState
{
	MGES_NORMAL = 0,
	MGES_CONSTRAINING = 1,
	MGES_CONSTRAINED = 2,
	MGES_FIXED = 3
};

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
	GeometricObject* 	m_pParent;
	VertexBase*			m_pVrtChild;
	byte				m_state;
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
	GeometricObject* 	m_pParent;
	VertexBase*			m_pVrtChild;
	EdgeBase* 			m_pEdgeChild[MG_EDGE_MAX_EDGE_CHILDREN];
	byte				m_numEdgeChildren;///< primarily required during refinement
	byte				m_state;
};

///	Holds information about face relations. Used internally.
struct MGFaceInfo
{
	MGFaceInfo()		{clear();}
	inline void clear()	{m_pParent = m_pVrtChild = NULL; m_numEdgeChildren = m_numFaceChildren = 0;}
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
	GeometricObject* 	m_pParent;
	VertexBase*			m_pVrtChild;
	EdgeBase* 			m_pEdgeChild[MG_FACE_MAX_EDGE_CHILDREN];
	Face*				m_pFaceChild[MG_FACE_MAX_FACE_CHILDREN];
	byte				m_numEdgeChildren;///< primarily required during refinement
	byte				m_numFaceChildren;///< primarily required during refinement
	byte				m_state;
};

///	Holds information about volume relations. Used internally.
struct MGVolumeInfo
{
	MGVolumeInfo()		{clear();}
	inline void clear()	{m_pParent = m_pVrtChild = NULL; m_numEdgeChildren = m_numFaceChildren = m_numVolChildren = 0;}
	inline bool has_children()	const {return m_pVrtChild || m_numEdgeChildren || m_numFaceChildren || m_numVolChildren;}
	inline void add_child(VertexBase* elem)	{assert(!m_pVrtChild); m_pVrtChild = elem;}
	inline void add_child(EdgeBase* elem)	{assert(m_numEdgeChildren < MG_VOLUME_MAX_EDGE_CHILDREN); m_pEdgeChild[m_numEdgeChildren++] = elem;}
	inline void add_child(Face* elem)		{assert(m_numFaceChildren < MG_VOLUME_MAX_FACE_CHILDREN); m_pFaceChild[m_numFaceChildren++] = elem;}
	inline void add_child(Volume* elem)		{assert(m_numVolChildren < MG_VOLUME_MAX_VOLUME_CHILDREN); m_pVolChild[m_numVolChildren++] = elem;}
	inline void remove_child(VertexBase* elem)	{m_pVrtChild = NULL;}
	inline void remove_child(EdgeBase* elem)	{m_numEdgeChildren = ArrayEraseEntry(m_pEdgeChild, elem, m_numEdgeChildren);}
	inline void remove_child(Face* elem)		{m_numFaceChildren = ArrayEraseEntry(m_pFaceChild, elem, m_numFaceChildren);}
	inline void remove_child(Volume* elem)		{m_numVolChildren = ArrayEraseEntry(m_pVolChild, elem, m_numVolChildren);}
	inline void replace_child(VertexBase* elem, VertexBase* child)	{assert(child == m_pVrtChild); m_pVrtChild = elem;}
	inline void replace_child(EdgeBase* elem, EdgeBase* child)		{ArrayReplaceEntry(m_pEdgeChild, elem, child, m_numEdgeChildren);}
	inline void replace_child(Face* elem, Face* child)				{ArrayReplaceEntry(m_pFaceChild, elem, child, m_numFaceChildren);}
	inline void replace_child(Volume* elem, Volume* child)			{ArrayReplaceEntry(m_pVolChild, elem, child, m_numVolChildren);}
	void unregister_from_children(MultiGrid& mg);

	GeometricObject* 	m_pParent;
	VertexBase*			m_pVrtChild;
	EdgeBase* 			m_pEdgeChild[MG_VOLUME_MAX_EDGE_CHILDREN];
	Face*				m_pFaceChild[MG_VOLUME_MAX_FACE_CHILDREN];
	Volume*				m_pVolChild[MG_VOLUME_MAX_VOLUME_CHILDREN];
	byte				m_numEdgeChildren;///< primarily required during refinement
	byte				m_numFaceChildren;///< primarily required during refinement
	byte				m_numVolChildren;///< primarily required during refinement
	byte				m_state;
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
 * In order to make state-assignment as effective as possible, one should
 * assign the status of the parent before creating its children.
 *
 * Note that states are currently not really used. They should probably
 * be removed.
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
		using Grid::get_geometric_object_collection;
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

		inline uint num_levels() const	{return m_hierarchy.num_subsets();}

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
		get_geometric_object_collection(int level)
		{return m_hierarchy.get_geometric_object_collection(level);}
		
	//	multi-level-geometric-object-collection
		virtual GeometricObjectCollection get_geometric_object_collection()
		{return m_hierarchy.get_geometric_object_collection();}
		
		template <class TElem> inline
		int get_level(TElem* elem)
		{return m_hierarchy.get_subset_index(elem);}

		GeometricObject* get_parent(GeometricObject* parent);

		template <class TElem> inline
		GeometricObject* get_parent(TElem* elem)
		{return get_info(elem).m_pParent;}

	//	number of children
		template <class TElem> inline
		bool has_children(TElem* elem) const
		{return get_info(elem).has_children();}

	////////////////////////////////
	//	CHILD QUANTITIES
	///	Returns the number of child vertices
		template <class TElem>
		inline size_t num_child_vertices(TElem* elem)	{return get_info(elem).m_pVrtChild ? 1 : 0;}

	///	Returns the number of child edges
	/**	\{	*/
		template <class TElem>
		inline size_t num_child_edges(TElem* elem)		{return get_info(elem).m_numEdgeChildren;}
		inline size_t num_child_edges(VertexBase*)		{return 0;}
	/**	\}	*/

	///	Returns the number of child faces
	/**	\{	*/
		template <class TElem>
		inline size_t num_child_faces(TElem* elem)		{return get_info(elem).m_numFaceChildren;}
		inline size_t num_child_faces(VertexBase*)		{return 0;}
		inline size_t num_child_faces(EdgeBase*)		{return 0;}
	/**	\}	*/

	///	Returns the number of child volumes
	/**	\{	*/
		inline size_t num_child_volumes(Volume* elem)		{return get_info(elem).m_numVolChildren;}
		template <class TElem>
		inline size_t num_child_volumes(TElem*)				{return 0;}
	/**	\}	*/


	////////////////////////////////
	//	CHILD ACCESS
	///	Returns the child vertex of the given element or NULL if there is none
		template <class TElem>
		inline VertexBase* get_child_vertex(TElem* elem)	{return get_info(elem).m_pVrtChild;}

	///	Returns the child edges of the given element or NULL if there is none
	/**	\{	*/
		template <class TElem>
		inline EdgeBase* get_child_edge(TElem* elem, size_t ind)	{return get_info(elem).m_pEdgeChild[ind];}
		inline EdgeBase* get_child_edge(VertexBase*, size_t)		{return NULL;}
	/**	\}	*/

	///	Returns the child faces of the given element or NULL if there is none
	/**	\{	*/
		template <class TElem>
		inline Face* get_child_face(TElem* elem, size_t ind)	{return get_info(elem).m_pFaceChild[ind];}
		inline Face* get_child_face(VertexBase*, size_t)		{return NULL;}
		inline Face* get_child_face(EdgeBase*, size_t)			{return NULL;}
	/**	\}	*/

	///	Returns the child volumes of the given element or NULL if there is none
	/**	\{	*/
		inline Volume* get_child_volume(Volume* elem, size_t ind)	{return get_info(elem).m_pVolChild[ind];}
		template <class TElem>
		inline Volume* get_child_volume(TElem*, size_t)	{return NULL;}
	/**	\}	*/

	////////////////////////////////
	//	STATUS ACCESS
	//	access to the elements multi-grid status
	///	returns one of the constants enumerated in MGElementStates.
		template <class TElem>
		inline byte get_status(TElem* elem)						{return get_info(elem).m_state;}
		
	///	changes the state of the vertex and adjusts the states of its children.
		void set_state(VertexBase* vrt, int state);
	///	changes the state of the edge and adjusts the states of its children.
		void set_state(EdgeBase* edge, int state);
	///	NOT YET IMPLEMENTED!
	/** changes the state of the face and adjusts the states of its children.*/
		void set_state(Face* face, int state);
	///	NOT YET IMPLEMENTED!
	/**	changes the state of the volume and adjusts the states of its children.*/
		void set_state(Volume* vol, int state);

	///	for debug purposes
		void check_edge_elem_infos(int level);
	///	for debug purposes
		void check_face_elem_infos(int level);
	///	for debug purposes
		void check_volume_elem_infos(int level);
		
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

		typedef Attachment<VertexInfo>	AVertexInfo;
		typedef Attachment<EdgeInfo>		AEdgeInfo;
		typedef Attachment<FaceInfo>		AFaceInfo;
		typedef Attachment<VolumeInfo>	AVolumeInfo;

	protected:
	//	initialization
		void init();
		
	//	info-access
		inline VertexInfo& get_info(VertexBase* v)	{return m_aaVrtInf[v];}
		inline EdgeInfo& get_info(EdgeBase* e)		{return m_aaEdgeInf[e];}
		inline FaceInfo& get_info(Face* f)			{return m_aaFaceInf[f];}
		inline VolumeInfo& get_info(Volume* v)		{return m_aaVolInf[v];}

	//	const info-access
		inline const VertexInfo& get_info(VertexBase* v) const	{return m_aaVrtInf[v];}
		inline const EdgeInfo& get_info(EdgeBase* e) const		{return m_aaEdgeInf[e];}
		inline const FaceInfo& get_info(Face* f) const			{return m_aaFaceInf[f];}
		inline const VolumeInfo& get_info(Volume* v) const		{return m_aaVolInf[v];}

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

	protected:
	//	hierarchy
		SubsetHandler	m_hierarchy;
		bool m_bHierarchicalInsertion;

	//	info attachments
		AVertexInfo	m_aVertexInfo;
		AEdgeInfo	m_aEdgeInfo;
		AFaceInfo	m_aFaceInfo;
		AVolumeInfo	m_aVolumeInfo;

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
