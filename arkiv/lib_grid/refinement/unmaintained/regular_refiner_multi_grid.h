/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_regular_refiner_multi_grid
#define __H__UG_regular_refiner_multi_grid

#include "lib_grid/tools/selector_multi_grid.h"
#include "refiner_interface.h"

namespace ug{

class UG_API RegularRefiner_MultiGrid : public IRefiner{
	public:
		RegularRefiner_MultiGrid ();
		RegularRefiner_MultiGrid (MultiGrid* pmg);

		void 		set_grid (MultiGrid* pmg);
		MultiGrid*	multi_grid () const;

		virtual Grid* get_associated_grid ();
		virtual Grid* grid ();

		virtual bool adaptivity_supported () const;
		virtual bool coarsening_supported () const;

		virtual bool mark (Vertex* v, RefinementMark refMark = RM_REFINE);
		virtual bool mark (Edge* e, RefinementMark refMark = RM_REFINE);
		virtual bool mark (Face* f, RefinementMark refMark = RM_REFINE);
		virtual bool mark (Volume* v, RefinementMark refMark = RM_REFINE);

		virtual void mark_neighborhood (size_t numIterations,
										RefinementMark refMark,
										bool sideNbrsOnly);

		virtual RefinementMark get_mark (Vertex* v);
		virtual RefinementMark get_mark (Edge* e);
		virtual RefinementMark get_mark (Face* f);
		virtual RefinementMark get_mark (Volume* v);

		virtual bool save_marks_to_file (const char* filename);

		template <typename TElem>
		bool is_closure (TElem* elem);

	protected:
		enum Marks{
			NONE = RM_NONE,
			COPY = RM_COPY,
			REGULAR = RM_REFINE,
			ANISOTROPIC = RM_ANISOTROPIC,
			LIFT = REGULAR | ANISOTROPIC | COPY
		};

		virtual void perform_refinement ();
		virtual bool perform_coarsening ();

	///	returns the number of locally marked edges on all levels of the hierarchy
		virtual void num_marked_edges_local (std::vector<int>& numMarkedEdgesOut);
	///	returns the number of locally marked faces on all levels of the hierarchy
		virtual void num_marked_faces_local (std::vector<int>& numMarkedFacesOut);
	///	returns the number of locally marked volumes on all levels of the hierarchy
		virtual void num_marked_volumes_local (std::vector<int>& numMarkedVolsOut);


		template <typename TElem>
		void adjust_side_states (
				size_t lvl,
				uint considerElemMarks,
				uint ignoreSideMarks,
				RefinementMark newSideMark,
				bool closure);

		template <typename TElem>
		void copy_state_to_sides (
				size_t lvl,
				uint considerElemMarks,
				bool closure);

		template <typename TSide>
		void adjust_side_of_states (
				size_t lvl,
				uint considerSideMarks,
				uint ignoreElemMarks,
				RefinementMark newElemMark,
				bool closure);

		template <typename TElem>
		void clear_dummies ();

		template <typename TElem>
		void mark_by_level_discrepancy (
				int lvl,
				Grid::VertexAttachmentAccessor<AInt> aaLvl);

		void collect_objects_for_refine ();

		bool refinement_is_allowed(Vertex* v);
		bool refinement_is_allowed(Edge* e);
		bool refinement_is_allowed(Face* f);
		bool refinement_is_allowed(Volume* v);

	private:
		RegularRefiner_MultiGrid (const RegularRefiner_MultiGrid&)	{}

		MGSelector						m_marks;
		MultiGrid*						m_pMG;

		ABool									m_aClosure;
		MultiElementAttachmentAccessor<AByte>	m_aaClosure;
};

}//	end of namespace

#endif