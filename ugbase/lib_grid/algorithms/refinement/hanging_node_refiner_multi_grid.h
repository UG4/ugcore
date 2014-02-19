// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 11.01.2011 (m,d,y)

#ifndef __H__UG__HANGIN_NODE_REFINER_MULTI_GRID__
#define __H__UG__HANGIN_NODE_REFINER_MULTI_GRID__

#include "hanging_node_refiner_base.h"
#include "lib_grid/tools/selector_multi_grid.h"

namespace ug
{

///	\addtogroup lib_grid_algorithms_refinement
///	@{

///	Specialization of ug::HangingNodeRefiner for ug::MultiGrid
/**	This class should be used, if hanging node refinement shall be
 * applied on a hierarchical grid (ug::MultiGrid).
 *
 * New elements will be constructed one level above their parent elements.
 *
 * HangingNodeRefiner_MultiGrid supports coarsening. Please note that coarsening
 * presumably has not yet reached its final implementation and behavior may
 * thus change in future revisions.
 * In the current implementation, refinement marks are removed during coarsening
 * and coarsening marks are removed during refinement. In order to use
 * coarsening and refinement on thus should first mark all elements which shall
 * be coarsened, perform coarsening, then mark all elements will shall be
 * refined and perform refinement (or vice versa).
 *
 * Take a look at ug::HangingNodeRefinerBase for a more in-depth documentation.
 *
 * \todo: 	Avoid the removal of refinement marks during coarsening and of
 * 			coarsening marks during refinement.
 *
 * \sa ug::HangingNodeRefinerBase, ug::HangingNodeRefiner_Grid
 */
class HangingNodeRefiner_MultiGrid : public HangingNodeRefinerBase<MGSelector>
{
	public:
		typedef HangingNodeRefinerBase<MGSelector>	BaseClass;
		using BaseClass::mark;

		enum HNodeCoarsenMarks{
			HNCM_FIRST = HNRM_MAX + 1,
			HNCM_NO_NBRS,
			HNCM_NONE,
			HNCM_PARTIAL,
			HNCM_REPLACE,
			HNCM_ALL,
			HNCM_INVALID,
			HNCM_UNKNOWN
		};

	public:
		HangingNodeRefiner_MultiGrid(IRefinementCallback* refCallback = NULL);
		HangingNodeRefiner_MultiGrid(MultiGrid& mg,
									IRefinementCallback* refCallback = NULL);

		virtual ~HangingNodeRefiner_MultiGrid();

		virtual void grid_to_be_destroyed(Grid* grid);

		virtual void assign_grid(MultiGrid& mg);
		virtual Grid* get_associated_grid()		{return m_pMG;}//depreciated
		virtual Grid* grid()					{return m_pMG;}
		virtual MultiGrid* multi_grid()			{return m_pMG;}

		virtual bool adaptivity_supported() const	{return true;}
		virtual bool coarsening_supported() const	{return true;}

	protected:
	///	performs coarsening on the elements marked with RM_COARSEN.
	/**
	 * The grid's message hub is informed using a "GridAdaption" message,
	 * passing an instance of GridMessage_Adapation, with values
	 * GMAT_HNODE_COARSENING_BEGINS and GMAT_HNODE_COARSENING_ENDS.
	 * See lib_grid/lib_grid_messages.h for more details.
	 *
	 * automatically adjusts the selection so that only valid coarsening operations
	 * are executed. Not that refinement marks are removed in this process.
	 *
	 * During a coarsening step only elements from the surface layer are removed.
	 * If not all children of a sub-surface element are marked, then the marks are ignored.
	 *
	 * Note that coarsen in contrary to refine is conservative. While refine
	 * extends the selection to construct a valid grid, coarsen shrinks the
	 * selection. On could think of implementing an alternative non conservative
	 * coarsen approach. However, the conservative one is the one mostly used
	 * in adaptive multigrid methods and was thus chosen here.
	 *
	 * coarsen returns false, if no elements have been coarsened, true if at
	 * least one has been coarsened.
	 */
		virtual bool perform_coarsening();

		void save_coarsen_marks_to_file(ISelector& sel, const char* filename);
		void debug_save(ISelector& sel, const char* filename);

	///	a callback that allows to deny refinement of special vertices
		virtual bool refinement_is_allowed(Vertex* elem);
	///	a callback that allows to deny refinement of special edges
		virtual bool refinement_is_allowed(EdgeBase* elem);
	///	a callback that allows to deny refinement of special faces
		virtual bool refinement_is_allowed(Face* elem);
	///	a callback that allows to deny refinement of special volumes
		virtual bool refinement_is_allowed(Volume* elem);

	///	performs registration and deregistration at a grid.
	/**	Initializes all grid related variables.
	 *  call set_grid(NULL) to unregister the observer from a grid.
	 *
	 * 	Please note that though the base grid features a set_grid method,
	 *  it is not declared virtual. This is because we want to call it
	 *  during construction and destruction.*/
		void set_grid(MultiGrid* mg);

	///	creates required vertices in higher levels.
		virtual void pre_refine();

	///	called before elements are removed in coarsening
	/**	Default implementation is emtpy */
		virtual void pre_coarsen()	{};

	///	called after elements have been removed in coarsening
	/**	Default implementation is emtpy */
		virtual void post_coarsen()	{};


	/**	Calls the base implementation and passes the mg-child vertices as
	 *  newCornerVrts. If newCornerVrts were passed to this method, they
	 *  are ignored.
	 *  \{*/
	///	calls base implementation and replaces cge with a normal edge.
		virtual void process_constraining_edge(ConstrainingEdge* cge);

		virtual void refine_edge_with_normal_vertex(EdgeBase* e,
											Vertex** newCornerVrts = NULL);
		virtual void refine_edge_with_hanging_vertex(EdgeBase* e,
											Vertex** newCornerVrts = NULL);

		virtual void refine_face_with_normal_vertex(Face* f,
											Vertex** newCornerVrts = NULL);
		virtual void refine_face_with_hanging_vertex(Face* f,
											Vertex** newCornerVrts = NULL);

		virtual void refine_volume_with_normal_vertex(Volume* v,
											Vertex** newVolumeVrts = NULL);
	/*	\} */

	///	Returns the vertex associated with the edge
		virtual Vertex* get_center_vertex(EdgeBase* e);

	///	Associates a vertex with the edge.
		virtual void set_center_vertex(EdgeBase* e, Vertex* v);

	///	Returns the vertex associated with the face
		virtual Vertex* get_center_vertex(Face* f);

	///	Associates a vertex with the face.
		virtual void set_center_vertex(Face* f, Vertex* v);


	///	collects corner vertices and fills them into the associated vector
	/**	The size of cornersOut is automatically adjusted.
	 *  The i-th element of corners out will contain the child vertex of the
	 *  i-th vertex of elem.*/
		template <class TElem>
		void collect_child_corners(std::vector<Vertex*>& cornersOut, TElem* elem)
		{
			cornersOut.resize(elem->num_vertices());
			for(size_t i = 0; i < elem->num_vertices(); ++i){
				//UG_ASSERT(m_pMG->get_child_vertex(elem->vertex(i)), "A child vertex has to exists!");
				cornersOut[i] = m_pMG->get_child_vertex(elem->vertex(i));
			}
		}

	///	Makes sure that only surface elements are marked and that only coarsen marks are used.
	/**	calls the overloaded template implementation for the four base element types.
	 * this method is called prior to restrict_selection_to_coarsen_families.*/
		virtual void restrict_selection_to_surface_coarsen_elements();

	///	Deselects all elements which have an unselected sibling.
	/** Siblings are elements who share the same parent.
	 * Calls the overloaded template implementation for the four base element types.
	 * this method is called after to restrict_selection_to_surface_coarsen_elements.*/
		virtual void restrict_selection_to_coarsen_families();

	///	Deselects all non-surface elements and all elements not marked with RM_COARSEN
		template <class TElem>
		void restrict_selection_to_surface_coarsen_elements();

	///	Only complete families (all siblings are selected) may be coarsened.
	/**	When calling this method, make sure that only surface elements are
	 * marked/selected and that only the mark RM_COARSEN is used.*/
		template <class TElem>
		void restrict_selection_to_coarsen_families();

	///	adjusts the selection marks to those specified in HNodeCoarsenMarks.
	/**	This method does not alter the selection itself, it only alters the
	 * selection-mark, associated with each selected element. It neither
	 * selects nor deselects any elements.
	 * It iterates over all selected elements of the given type and assigns one
	 * of the marks, depending on how many associated elements of the next higher
	 * dimension are selected.
	 * Note that this method should only be called for vertices, edges and faces.
	 * It does not make sense to call it for volumes.*/
//		template <class TElem>
//		void classify_selection();

	///	Applies marks like RM_COARSEN_CONSTRAINING or RM_COARSEN_UNCONSTRAIN
	/**	This method should first be called for Face, then for EdgeBase, then for
	 * Vertex.*/
		//template <class TElem>
		//void adjust_coarsen_marks_on_side_elements();

	///deselect coarsen families, which are adjacent to unselected constraining elements
	/**	If at least one family was deselected, this method returns true.*/
//		template <class TElem>
//		bool deselect_invalid_coarsen_families();
//
//		template <class TElem>
//		void deselect_isolated_sides();
//
//		template <class TElem>
//		void deselect_uncoarsenable_parents();

	///	called by the coarsen method in order to adjust the selection to valid elements.
	/**	This method is responsible to mark all elements that shall be coarsened.
	 * Only sub-surface elements may be coarsened. If a sub-surface element has
	 * a side, which is not a sub-surface element, then the element may not be
	 * coarsened.
	 * If scheduleCoarseningBeginsMessage is set to true, the method will schedule
	 * a GridMessage_Adaption message which indicates that coarsening begins.*/
		virtual void collect_objects_for_coarsen(bool scheduleCoarseningBeginsMessage = false);


	////////////////////////////////////////
	//	Callbacks for PARALLELIZATION

	///	called to check, whether another iteration of collect_objects_for_coarsen has to be performed.
	/**	This method is especially useful if coarsening shall be applied in a
	 * parallel environment. Each process calls this method and basses a boolean,
	 * which indicates if it requires the continuation of collect_objects_for_coarsen.
	 * If one process does require continuation, then all processes have to continue.
	 * The default implementation simply returns the local status.*/
//		virtual bool continue_collect_objects_for_coarsen(bool continueRequired)
//		{return continueRequired;}

	///	called during collect_objects_for_coarsen when coarsen-marks shall be distributed
	/**	This method is especially useful if coarsening shall be applied in a
	 * parallel environment. A derived class can implement this method to schedule
	 * communication. The default implementation does nothing.
	 * \{ */
//		virtual void broadcast_vertex_coarsen_marks()	{}
//		virtual void broadcast_edge_coarsen_marks()		{}
//		virtual void broadcast_face_coarsen_marks()		{}
	/** \} */

	///	allows to check whether a distributed grid contains edges
	/**	The default implementation returns whether the local grid contains edges.*/
		virtual bool contains_edges()			{return m_pMG->num<EdgeBase>() > 0;}

	///	allows to check whether a distributed grid contains faces
	/**	The default implementation returns whether the local grid contains faces.*/
		virtual bool contains_faces()			{return m_pMG->num<Face>() > 0;}

	///	allows to check whether a distributed grid contains volumes
	/**	The default implementation returns whether the local grid contains volumes.*/
		virtual bool contains_volumes()			{return m_pMG->num<Volume>() > 0;}

	/**	This callback is called during execution of the coarsen() method after
	 * collect_objects_for_coarsen is done. It is responsible to mark
	 * elements for hnode coarsening. That means all elements on which a hanging
	 * node or constrained children shall be created have to be marked using
	 * mark_for_hnode_refinement during this method. The default implementation
	 * performs this marking for all local elements.*/
		//virtual void assign_hnode_coarsen_marks();

		virtual void broadcast_marks_horizontally(bool vertices, bool edges, bool faces,
												  bool allowDeselection = false)	{}
		virtual void broadcast_marks_vertically(bool vertices, bool edges,
												bool faces, bool volumes,
												bool allowDeselection = false)	{}

		virtual void copy_marks_to_vmasters(bool vertices, bool edges,
											bool faces, bool volumes)			{}

		virtual void copy_marks_to_vslaves(bool vertices, bool edges,
											bool faces, bool volumes)			{}

		virtual bool one_proc_true(bool localProcTrue)		{return localProcTrue;}

	private:
		MultiGrid*	m_pMG;
};

/// @}	// end of add_to_group command

}//	end of namespace

#endif
