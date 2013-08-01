// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 11.01.2011 (m,d,y)

#ifndef __H__UG__ADAPTIVE_REGULAR_REFINER_MULTI_GRID__
#define __H__UG__ADAPTIVE_REGULAR_REFINER_MULTI_GRID__

#include "hanging_node_refiner_multi_grid.h"
#include "lib_grid/tools/selector_grid.h"

namespace ug
{

///	\addtogroup lib_grid_algorithms_refinement
///	@{

///	Specialization of IRefiner for adaptive multigrid refinement with closure.
/**	The implementation is based on HangingNodeRefiner_MultiGrid.
 * Hanging nodes are most likely present in the resulting grids, even though
 * they should not be used by discretizations etc. Their sole purpose is to
 * allow efficient grid adaption
 *
 * \todo	The current implementation has a performance overhead, since currently
 * 			all closure elements are removed prior to grid adaption. This should be
 * 			reduced to only remove closure elements which are affected by adaption.
 * 			One could therefore replace Selector by BoolMarker and adjust the
 * 			collect_objects_for_refinement / coarsening routines to use callbacks
 * 			when identifying valid refine / coarsen objects.
 */
class AdaptiveRegularRefiner_MultiGrid : public HangingNodeRefiner_MultiGrid
{
	public:
		typedef HangingNodeRefiner_MultiGrid BaseClass;
		using HangingNodeRefiner_MultiGrid::mark;

	public:
		AdaptiveRegularRefiner_MultiGrid(IRefinementCallback* refCallback = NULL);
		AdaptiveRegularRefiner_MultiGrid(MultiGrid& mg,
										 IRefinementCallback* refCallback = NULL);

		virtual ~AdaptiveRegularRefiner_MultiGrid();

		virtual void assign_grid(MultiGrid& mg);

	protected:
	///	performs registration and deregistration at a grid.
	/**	Initializes all grid related variables.
	 *  call set_grid(NULL) to unregister the observer from a grid.
	 *
	 * 	Please note that though the base grid features a set_grid method,
	 *  it is not declared virtual. This is because we want to call it
	 *  during construction and destruction.*/
		void set_grid(MultiGrid* mg);

	///	remove closure elements
		void remove_closure_elements();

	///	calls create_closure_elements_2d / 3d, depending on the presence of volume elements.
		void create_closure_elements();

	///	creates closure elements for 2d geometries
		void create_closure_elements_2d();

	///	creates closure elements for 3d geometries
		void create_closure_elements_3d();

	///	collects parents of all closure elements which share a value with the given mark.
	/**	Note that the parents container is not cleared. The method may thus be called
	 * multiple times.
	 * If a parent is a closure element (which shouldn't be the case) it won't be
	 * added to the container.
	 * Note that a parent may be added multiple times to the parents container.*/
		template <class TElem>
		void get_parents_of_marked_closure_elements(std::vector<GeometricObject*>& parents,
									   	   	   	    Selector::status_t mark);

	///	removes all closure elements, calls the base implementation and creates a new closure
		virtual void perform_refinement();

	///	removes all closure elements, calls the base implementation and creates a new closure
		virtual bool perform_coarsening();

	protected:
		Selector	m_closureElems;
};

/// @}	// end of add_to_group command

}//	end of namespace

#endif
