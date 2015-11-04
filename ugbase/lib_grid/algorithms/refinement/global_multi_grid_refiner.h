#ifndef __H__LIB_GRID__GLOBAL_MULTI_GRID_REFINER__
#define __H__LIB_GRID__GLOBAL_MULTI_GRID_REFINER__

#include <vector>
#include "lib_grid/lg_base.h"
#include "lib_grid/multi_grid.h"
#include "refinement_callbacks.h"
#include "refiner_interface.h"

namespace ug
{

///	\addtogroup lib_grid_algorithms_refinement
///	@{

class GlobalMultiGridRefiner : public IRefiner, public GridObserver
{
	public:
		GlobalMultiGridRefiner(IRefinementCallback* refCallback = NULL);
		GlobalMultiGridRefiner(MultiGrid& mg,
							   IRefinementCallback* refCallback = NULL);
							   
		virtual ~GlobalMultiGridRefiner();

		virtual void grid_to_be_destroyed(Grid* grid);
		
		void assign_grid(MultiGrid& mg);
		void assign_grid(MultiGrid* mg);

		virtual Grid* get_associated_grid()		{return m_pMG;}
		virtual Grid* grid()					{return m_pMG;}
		virtual MultiGrid* multi_grid()			{return m_pMG;}

		virtual bool adaptivity_supported() const	{return false;}
		virtual bool coarsening_supported() const	{return false;}

		virtual bool save_marks_to_file(const char* filename);

	protected:
	///	returns the number of (globally) marked edges on this level of the hierarchy
		virtual void num_marked_edges_local(std::vector<int>& numMarkedEdgesOut);
	///	returns the number of (globally) marked faces on this level of the hierarchy
		virtual void num_marked_faces_local(std::vector<int>& numMarkedFacesOut);
	///	returns the number of (globally) marked volumes on this level of the hierarchy
		virtual void num_marked_volumes_local(std::vector<int>& numMarkedVolsOut);

		template <class TElem>
		void num_marked_elems(std::vector<int>& numMarkedElemsOut);

	////////////////////////////////
	///	performs refinement on the marked elements.
		virtual void perform_refinement();

	///	a callback that allows to deny refinement of special vertices
		virtual bool refinement_is_allowed(Vertex* elem)	{return true;}
	///	a callback that allows to deny refinement of special edges
		virtual bool refinement_is_allowed(Edge* elem)		{return true;}
	///	a callback that allows to deny refinement of special faces
		virtual bool refinement_is_allowed(Face* elem)			{return true;}
	///	a callback that allows to deny refinement of special volumes
		virtual bool refinement_is_allowed(Volume* elem)		{return true;}
		
	///	this method helps derived classes to perform operations directly before actual element refinment is performed.
	/**	Called from the refine() method in each refinement-iteration after
	 *	collect_objects_for_refine().
	 *	Default implementation is empty.*/
		virtual void refinement_step_begins()	{};

	///	this method helps derived classes to perform operations directly after actual element refinment took place.
	/**	Called from the refine() method in each refinement-iteration after
	 *	all scheduled elements had been refined.
	 *	The refine process will either terminate after this method or will
	 *	start a new iteration, if new elements had been marked during refine.
	 *	Default implementation is empty.*/
		virtual void refinement_step_ends()		{};
		
	protected:
		MultiGrid*	m_pMG;
};

/// @}
}//	end of namespace

#endif
