/*
 * Copyright (c) 2014-2019:  G-CSC, Goethe University Frankfurt
 * Author: Martin Stepniewski
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

#ifndef __H__LIB_GRID__GLOBAL_SUBDIVISION_MULTI_GRID_REFINER__
#define __H__LIB_GRID__GLOBAL_SUBDIVISION_MULTI_GRID_REFINER__

#include "lib_grid/lg_base.h"
#include "lib_grid/multi_grid.h"
#include "global_multi_grid_refiner.h"
#include "lib_grid/algorithms/subdivision/subdivision_volumes.h"

namespace ug
{

///	\addtogroup lib_grid_algorithms_refinement
///	@{

/// Specialization of the GlobalMultiGridRefiner class to incorporate subdivision refinement
/**	Each step of global refinement is subsequently followed by
 * a vertex smoothing pass determined by a user-specified subdivision scheme.
 *
 * Currently this only works in 2d and 3d.
 */
template <typename TAPosition>
class GlobalSubdivisionMultiGridRefiner : public GlobalMultiGridRefiner
{
	public:
		GlobalSubdivisionMultiGridRefiner(SPRefinementProjector projector = nullptr);

		GlobalSubdivisionMultiGridRefiner(MultiGrid& mg, SPRefinementProjector projector = nullptr);

		GlobalSubdivisionMultiGridRefiner(MultiGrid& mg, TAPosition& aPos, MGSubsetHandler& sh,
								MGSubsetHandler& markSH, SPRefinementProjector projector = nullptr);

		~GlobalSubdivisionMultiGridRefiner() override;

	///	sets the manifold subsets which shall be linearly refined
		void set_linear_manifold_subsets(MGSubsetHandler& linearManifoldSH, const char* linearManifoldSubsets);

	///	sets the default SubsetHandler
		void assign_subset_handler(MGSubsetHandler* sh) {m_pSH = sh;}
		void assign_subset_handler(MGSubsetHandler& sh) {m_pSH = &sh;}

	///	sets the SubsetHandler designated for obligatory marked manifold elements
		void assign_mark_subset_handler(MGSubsetHandler* markSH) {m_pMarkSH = markSH;}
		void assign_mark_subset_handler(MGSubsetHandler& markSH) {m_pMarkSH = &markSH;}

	///	sets the position attachment depending on the world dimension
		void assign_position_attachment(TAPosition* aPos);
		void assign_position_attachment(TAPosition& aPos);

	///	sets constrained subdivision volumes scheme
		void set_constrained_subdivision(bool constrained) {m_bConstrained = constrained;}

	///	projection of the vertices of all levels to their smooth subdivision limit positions to ensure node nested hierarchy
		void nest_hierarchy();

	protected:
	///	performs subdivision smoothing on the marked elements after base class regular refinement
		void smooth();

	///	wrapper for smooth() method in parallel case (see class ParallelGlobalSubdivisionRefiner)
		void refinement_step_ends() override;

	protected:
		MGSubsetHandler* m_pSH;
		MGSubsetHandler* m_pMarkSH;
		MGSubsetHandler* m_spLinearManifoldSH;
//		SmartPtr<MGSubsetHandler> m_spLinearManifoldSH;
		TAPosition* m_pAPos;
		bool m_bConstrained;
};

/// @}
}//	end of namespace

#endif
