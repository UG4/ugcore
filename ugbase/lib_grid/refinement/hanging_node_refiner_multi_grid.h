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
		using BaseClass = HangingNodeRefinerBase;
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
		explicit HangingNodeRefiner_MultiGrid(SPRefinementProjector projector = nullptr);

		explicit HangingNodeRefiner_MultiGrid(MultiGrid& mg, SPRefinementProjector projector = nullptr);

		~HangingNodeRefiner_MultiGrid() override;

		void grid_to_be_destroyed(Grid* grid) override;

		virtual void assign_grid(MultiGrid& mg);

		Grid* get_associated_grid() override {return m_pMG;}//depreciated
		Grid* grid() override {return m_pMG;}
		virtual MultiGrid* multi_grid() {return m_pMG;}

		bool adaptivity_supported() const override {return true;}
		bool coarsening_supported() const override {return true;}

	protected:
	///	returns the number of (globally) marked edges on this level of the hierarchy
		void num_marked_edges_local(std::vector<int>& numMarkedEdgesOut) override;
	///	returns the number of (globally) marked faces on this level of the hierarchy
		void num_marked_faces_local(std::vector<int>& numMarkedFacesOut) override;
	///	returns the number of (globally) marked volumes on this level of the hierarchy
		void num_marked_volumes_local(std::vector<int>& numMarkedVolsOut) override;
		
		template <typename TElem>
		void num_marked_elems(std::vector<int>& numMarkedElemsOut);
		
	///	performs coarsening on the elements marked with RM_COARSEN.
	/**
	 * The grid's message hub is informed using a "GridAdaption" message,
	 * passing an instance of GridMessage_Adaption, with values
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
	 * selection. On could think of implementing an alternative non-conservative
	 * coarsen approach. However, the conservative one is the one mostly used
	 * in adaptive multigrid methods and was thus chosen here.
	 *
	 * coarsen returns false, if no elements have been coarsened, true if at
	 * least one has been coarsened.
	 */
		bool perform_coarsening() override;

		void save_coarsen_marks_to_file(ISelector& sel, const char* filename);
		void debug_save(ISelector& sel, const char* filename);

	///	a callback that allows to deny refinement of special vertices
		bool refinement_is_allowed(Vertex* elem) override;
	///	a callback that allows to deny refinement of special edges
		bool refinement_is_allowed(Edge* elem) override;
	///	a callback that allows to deny refinement of special faces
		bool refinement_is_allowed(Face* elem) override;
	///	a callback that allows to deny refinement of special volumes
		bool refinement_is_allowed(Volume* elem) override;

	///	performs registration and deregistration at a grid.
	/**	Initializes all grid related variables.
	 *  call set_grid(nullptr) to unregister the observer from a grid.
	 *
	 * 	Please note that though the base grid features a set_grid method,
	 *  it is not declared virtual. This is because we want to call it
	 *  during construction and destruction.*/
		void set_grid(MultiGrid* mg);

		void assign_hnode_marks() override;

	///	creates required vertices in higher levels.
		void pre_refine() override;

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
		void process_constraining_edge(ConstrainingEdge* cge) override;

		void refine_edge_with_normal_vertex(Edge* e, Vertex** newCornerVrts = nullptr) override;

		void refine_edge_with_hanging_vertex(Edge* e, Vertex** newCornerVrts = nullptr) override;

		void refine_face_with_normal_vertex(Face* f, Vertex** newCornerVrts = nullptr) override;

		void refine_face_with_hanging_vertex(Face* f, Vertex** newCornerVrts = nullptr) override;

		void refine_volume_with_normal_vertex(Volume* v, Vertex** newVolumeVrts = nullptr) override;
	/*	\} */

	///	Returns the vertex associated with the edge
		Vertex* get_center_vertex(Edge* e) override;

	///	Associates a vertex with the edge.
		void set_center_vertex(Edge* e, Vertex* v) override;

	///	Returns the vertex associated with the face
		Vertex* get_center_vertex(Face* f) override;

	///	Associates a vertex with the face.
		void set_center_vertex(Face* f, Vertex* v) override;


	///	collects corner vertices and fills them into the associated vector
	/**	The size of cornersOut is automatically adjusted.
	 *  The i-th element of corners out will contain the child vertex of the
	 *  i-th vertex of elem.*/
		template <typename TElem>
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
		template <typename TElem>
		void restrict_selection_to_surface_coarsen_elements();

	///	Only complete families (all siblings are selected) may be coarsened.
	/**	When calling this method, make sure that only surface elements are
	 * marked/selected and that only the mark RM_COARSEN is used.*/
		template <typename TElem>
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
//		template <typename TElem>
//		void classify_selection();

	///	Applies marks like RM_COARSEN_CONSTRAINING or RM_COARSEN_UNCONSTRAIN
	/**	This method should first be called for Face, then for Edge, then for
	 * Vertex.*/
		//template <typename TElem>
		//void adjust_coarsen_marks_on_side_elements();

	///deselect coarsen families, which are adjacent to unselected constraining elements
	/**	If at least one family was deselected, this method returns true.*/
//		template <typename TElem>
//		bool deselect_invalid_coarsen_families();
//
//		template <typename TElem>
//		void deselect_isolated_sides();
//
//		template <typename TElem>
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
		virtual bool contains_edges() {return m_pMG->num<Edge>() > 0;}

	///	allows to check whether a distributed grid contains faces
	/**	The default implementation returns whether the local grid contains faces.*/
		virtual bool contains_faces() {return m_pMG->num<Face>() > 0;}

	///	allows to check whether a distributed grid contains volumes
	/**	The default implementation returns whether the local grid contains volumes.*/
		virtual bool contains_volumes() {return m_pMG->num<Volume>() > 0;}

	/**	This callback is called during execution of the coarsen() method after
	 * collect_objects_for_coarsen is done. It is responsible to mark
	 * elements for hnode coarsening. That means all elements on which a hanging
	 * node or constrained children shall be created have to be marked using
	 * mark_for_hnode_refinement during this method. The default implementation
	 * performs this marking for all local elements.*/
		//virtual void assign_hnode_coarsen_marks();

		virtual void broadcast_marks_horizontally(bool vertices, bool edges, bool faces,
												  bool allowDeselection = false) {}

		virtual void broadcast_marks_vertically(bool vertices, bool edges,
												bool faces, bool volumes,
												bool allowDeselection = false) {}

		virtual void copy_marks_to_vmasters(bool vertices, bool edges, bool faces, bool volumes) {}

		virtual void copy_marks_to_vslaves(bool vertices, bool edges, bool faces, bool volumes) {}

		virtual bool one_proc_true(bool localProcTrue) {return localProcTrue;}

	private:
		MultiGrid*	m_pMG;
};

/// @}	// end of add_to_group command

}//	end of namespace

#endif
