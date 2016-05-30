/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_GRID__GLOBAL_MULTI_GRID_REFINER__
#define __H__LIB_GRID__GLOBAL_MULTI_GRID_REFINER__

#include <vector>
#include "lib_grid/lg_base.h"
#include "lib_grid/multi_grid.h"
#include "refiner_interface.h"

namespace ug
{

///	\addtogroup lib_grid_algorithms_refinement
///	@{

class GlobalMultiGridRefiner : public IRefiner, public GridObserver
{
	public:
		GlobalMultiGridRefiner(SPRefinementProjector projector = SPNULL);
		GlobalMultiGridRefiner(MultiGrid& mg,
							   SPRefinementProjector projector = SPNULL);
							   
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
