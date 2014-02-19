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
#include "multi_grid_child_info.h"
#include "algorithms/attachment_util.h"

namespace ug
{

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
		using Grid::get_grid_objects;
		using Grid::create;
		using Grid::create_by_cloning;


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
	 * TGeomObj has to be a geometric object type as described in grid_base_objects.h.
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
		VertexIterator create_by_cloning(Vertex* pCloneMe,
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
		inline GridObjectCollection
		get_grid_objects(int level)
		{return m_hierarchy.get_grid_objects_in_subset(level);}
		
	//	multi-level-geometric-object-collection
		virtual GridObjectCollection get_grid_objects()
		{return m_hierarchy.get_grid_objects();}
		
		template <class TElem> inline
		int get_level(TElem* elem) const
		{return m_hierarchy.get_subset_index(elem);}

		GridObject* get_parent(GridObject* parent) const;
		inline GridObject* get_parent(Vertex* o) const	{return get_info(o).m_pParent;}
		inline GridObject* get_parent(EdgeBase* o) const	{return get_info(o).m_pParent;}
		inline GridObject* get_parent(Face* o) const		{return m_aaParentFACE[o];}
		inline GridObject* get_parent(Volume* o) const		{return m_aaParentVOL[o];}

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
		size_t num_children(GridObject* elem)	const;
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
		inline size_t num_child_edges(Vertex*) const	{return 0;}
	/**	\}	*/

	///	Returns the number of child faces
	/**	\{	*/
		template <class TElem>
		inline size_t num_child_faces(TElem* elem) const	{return get_info(elem).num_child_faces();}
		inline size_t num_child_faces(Vertex*) const	{return 0;}
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
		TChild* get_child(GridObject* elem, size_t ind) const;
	/** \} */

	///	Returns the child vertex of the given element or NULL if there is none
		template <class TElem>
		inline Vertex* get_child_vertex(TElem* elem) const	{return get_info(elem).child_vertex();}

	///	Returns the child edges of the given element or NULL if there is none
	/**	\{	*/
		template <class TElem>
		inline EdgeBase* get_child_edge(TElem* elem, size_t ind) const	{return get_info(elem).child_edge(ind);}
		inline EdgeBase* get_child_edge(Vertex*, size_t) const		{return NULL;}
	/**	\}	*/

	///	Returns the child faces of the given element or NULL if there is none
	/**	\{	*/
		template <class TElem>
		inline Face* get_child_face(TElem* elem, size_t ind) const	{return get_info(elem).child_face(ind);}
		inline Face* get_child_face(Vertex*, size_t) const		{return NULL;}
		inline Face* get_child_face(EdgeBase*, size_t) const			{return NULL;}
	/**	\}	*/

	///	Returns the child volumes of the given element or NULL if there is none
	/**	\{	*/
		inline Volume* get_child_volume(Volume* elem, size_t ind) const	{return get_info(elem).child_volume(ind);}
		template <class TElem>
		inline Volume* get_child_volume(TElem*, size_t) const	{return NULL;}
	/**	\}	*/

	///	clears the relation between a parent and its children
	/**	Use with care. This method should only be called if no other option exists.*/
		template <class TElem>
		void clear_child_connections(TElem* parent);

	///	establishes a parent child connection between the given elements
	/**	Note that the MultGrid class automatically establishes a parent child
	 * connection during element creation, based on the suppliend parent argument.
	 * This method should thus only be used in the rare cases where this
	 * automatic association is not sufficient.
	 * Note that only elements of equal or higher dimension can be parent to a
	 * given element.
	 * Note that while parent may be NULL, elem has to be supplied.
	 * The method also sets the parent type if a parent is supplied and leaves the
	 * parent type as it is if none is supplied.*/
		template <class TElem>
		void associate_parent(TElem* elem, GridObject* parent);

	///	returns the object-type of the parent of a given object
		template <class TElem>
		char parent_type(TElem* elem) const;

	///	sets the object-type of the parent of a given object
	/**	The parent type is normally handled internally. However, e.g. during
	 * parallel redistribution it may have to be set from outside (e.g. if
	 * a parent element hasn't been transfered to the same process as its children).*/
		template <class TElem>
		void set_parent_type(TElem* elem, char type);

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
		virtual void vertex_created(Grid* grid, Vertex* vrt,
									GridObject* pParent = NULL,
									bool replacesParent = false);

	 /**  In order to correctly register e in the hierarchy, we have to
	 *  replace pParent with e in the list of children of pParents parent.
	 *  This means that if a grid-observer registered after the multi-grid itself,
	 *  iterates over the list of children of pParents parent in its
	 *  edge_created method, it won't find pParent.*/
		virtual void edge_created(Grid* grid, EdgeBase* e,
									GridObject* pParent = NULL,
									bool replacesParent = false);

	 /**  In order to correctly register f in the hierarchy, we have to
	 *  replace pParent with f in the list of children of pParents parent.
	 *  This means that if a grid-observer registered after the multi-grid itself,
	 *  iterates over the list of children of pParents parent in its
	 *  face_created method, it won't find pParent.*/
		virtual void face_created(Grid* grid, Face* f,
									GridObject* pParent = NULL,
									bool replacesParent = false);

	 /**  In order to correctly register vol in the hierarchy, we have to
	 *  replace pParent with vol in the list of children of pParents parent.
	 *  This means that if a grid-observer registered after the multi-grid itself,
	 *  iterates over the list of children of pParents parent in its
	 *  volume_created method, it won't find pParent.*/
		virtual void volume_created(Grid* grid, Volume* vol,
									GridObject* pParent = NULL,
									bool replacesParent = false);

		virtual void vertex_to_be_erased(Grid* grid, Vertex* vrt,
										 Vertex* replacedBy = NULL);

		virtual void edge_to_be_erased(Grid* grid, EdgeBase* e,
										 EdgeBase* replacedBy = NULL);

		virtual void face_to_be_erased(Grid* grid, Face* f,
										 Face* replacedBy = NULL);

		virtual void volume_to_be_erased(Grid* grid, Volume* vol,
										 Volume* replacedBy = NULL);

	protected:

	//	Note: VertexInfo and EdgeInfo are stored directly, FaceInfo and
	//	VolumeInfo are stored dynamically.
		typedef Attachment<GridObject*>	AParent;
		typedef Attachment<VertexInfo>			AVertexInfo;
		typedef Attachment<EdgeInfo>			AEdgeInfo;
		typedef Attachment<FaceInfo*>			AFaceInfo;
		typedef Attachment<VolumeInfo*>			AVolumeInfo;
		typedef Attachment<char>				AParentType;


	//	initialization
		void init();
		
	//	create levels
		void create_levels(int numLevels);

	//	info-access
		inline VertexInfo& get_info(Vertex* v);
		inline EdgeInfo& get_info(EdgeBase* e);
		inline FaceInfo& get_info(Face* f);
		inline VolumeInfo& get_info(Volume* v);

	//	const info-access
		inline const VertexInfo& get_info(Vertex* v) const;
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
		inline size_t num_children(TElem* elem, const Vertex&) const
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
		inline Vertex* get_child(TElem* elem, size_t ind, const Vertex&) const
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
		inline void set_parent(Vertex* o, GridObject* p)	{get_info(o).m_pParent = p;}
		inline void set_parent(EdgeBase* o, GridObject* p)		{get_info(o).m_pParent = p;}
		inline void set_parent(Face* o, GridObject* p)			{m_aaParentFACE[o] = p;}
		inline void set_parent(Volume* o, GridObject* p)		{m_aaParentVOL[o] = p;}
	/**	\} */

	///	adds a child to the given object
	/** \{ */
		template <class TParent, class TChild>
		void add_child(TParent* p, TChild* c);

		template <class TChild>
		void add_child(GridObject* p, TChild* c);
	/** \} */

	///	removes a child from the given object
	/** \{ */
		template <class TParent, class TChild>
		void remove_child(TParent* p, TChild* c);

		template <class TChild>
		void remove_child(GridObject* p, TChild* c);
	/** \} */

	///	creates the info-object for the given object (if necessary)
	/**	\{ */
		inline void create_child_info(Vertex* o){}
		inline void create_child_info(EdgeBase* o)	{}
		inline void create_child_info(Face* o)		{if(!m_aaFaceInf[o]) m_aaFaceInf[o] = new FaceInfo();}
		inline void create_child_info(Volume* o)	{if(!m_aaVolInf[o]) m_aaVolInf[o] = new VolumeInfo();}
	/**	\} */

	///	releases the info-object for the given object (if necessary)
	/**	\{ */
		inline void release_child_info(Vertex* o)	{}
		inline void release_child_info(EdgeBase* o)		{}
		inline void release_child_info(Face* o)			{if(m_aaFaceInf[o]) delete m_aaFaceInf[o]; m_aaFaceInf[o] = NULL;}
		inline void release_child_info(Volume* o)		{if(m_aaVolInf[o]) delete m_aaVolInf[o]; m_aaVolInf[o] = NULL;}
	/**	\} */


	//	hierarchy
		SubsetHandler	m_hierarchy;
		bool m_bHierarchicalInsertion;

	//	parent attachment
		AParent		m_aParent;

	//	info attachments
		AVertexInfo	m_aVertexInfo;
		AEdgeInfo	m_aEdgeInfo;
		AFaceInfo	m_aFaceInfo;
		AVolumeInfo	m_aVolumeInfo;
		AParentType	m_aParentType;

	//	parent access - only required for faces and volumes.
		Grid::FaceAttachmentAccessor<AParent>		m_aaParentFACE;
		Grid::VolumeAttachmentAccessor<AParent>		m_aaParentVOL;

	//	element info access
		Grid::VertexAttachmentAccessor<AVertexInfo>	m_aaVrtInf;
		Grid::EdgeAttachmentAccessor<AEdgeInfo>		m_aaEdgeInf;
		Grid::FaceAttachmentAccessor<AFaceInfo>		m_aaFaceInf;
		Grid::VolumeAttachmentAccessor<AVolumeInfo>	m_aaVolInf;

		MultiElementAttachmentAccessor<AParentType>	m_aaParentType;
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
