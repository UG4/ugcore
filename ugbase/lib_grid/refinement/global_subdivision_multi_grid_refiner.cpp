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

#include "global_subdivision_multi_grid_refiner.h"


using namespace std;

namespace ug
{

template <class TAPosition>
GlobalSubdivisionMultiGridRefiner<TAPosition>::
GlobalSubdivisionMultiGridRefiner(SPRefinementProjector projector) :
	GlobalMultiGridRefiner(projector)
{
	m_pMG = NULL;
	m_pSH = NULL;
	m_pMarkSH = NULL;
	m_spLinearManifoldSH = NULL;
	m_pAPos = NULL;
	m_bConstrained = false;
}

template <class TAPosition>
GlobalSubdivisionMultiGridRefiner<TAPosition>::
GlobalSubdivisionMultiGridRefiner(MultiGrid& mg, SPRefinementProjector projector) :
	GlobalMultiGridRefiner(mg, projector)
{
	m_pSH = NULL;
	m_pMarkSH = NULL;
	m_spLinearManifoldSH = NULL;
	m_pAPos = NULL;
	m_bConstrained = false;
}

template <class TAPosition>
GlobalSubdivisionMultiGridRefiner<TAPosition>::
GlobalSubdivisionMultiGridRefiner(MultiGrid& mg, TAPosition& aPos, MGSubsetHandler& sh,
							MGSubsetHandler& markSH,
							SPRefinementProjector projector) :
	GlobalMultiGridRefiner(mg, projector)
{
	m_pSH = NULL;
	assign_subset_handler(sh);

	m_pMarkSH = NULL;
	assign_mark_subset_handler(markSH);

//	m_spLinearManifoldSH = SmartPtr<MGSubsetHandler>(new MGSubsetHandler(*m_pMG));
	m_spLinearManifoldSH = NULL;
//	set_linear_manifold_subsets(*m_spLinearManifoldSH, "");

	m_pAPos = NULL;
	assign_position_attachment(aPos);
}

template <class TAPosition>
GlobalSubdivisionMultiGridRefiner<TAPosition>::~GlobalSubdivisionMultiGridRefiner()
{
	if(m_pMG)
		m_pMG->unregister_observer(this);
}

template <class TAPosition>
void GlobalSubdivisionMultiGridRefiner<TAPosition>::set_linear_manifold_subsets(MGSubsetHandler& linearManifoldSH,
																				const char* linearManifoldSubsets)
{
	m_spLinearManifoldSH = &linearManifoldSH;
	InitLinearManifoldSubsetHandler(*m_pMG, *m_pSH, *m_spLinearManifoldSH, linearManifoldSubsets);
}

template <class TAPosition>
void GlobalSubdivisionMultiGridRefiner<TAPosition>::assign_position_attachment(TAPosition& aPos)
{
	if(TAPosition::ValueType::Size == 1){
		UG_THROW("Error in GlobalSubdivisionMultiGridRefiner::assign_position_attachment:\n"
				 "Currently only dimensions 2 and 3 are supported.\n");
	}
	else
		m_pAPos = &aPos;
}

template <class TAPosition>
void GlobalSubdivisionMultiGridRefiner<TAPosition>::assign_position_attachment(TAPosition* aPos)
{
	if(TAPosition::ValueType::Size == 1){
		UG_THROW("Error in GlobalSubdivisionMultiGridRefiner::assign_position_attachment:\n"
				 "Currently only dimensions 2 and 3 are supported.\n");
	}
	else
		m_pAPos = aPos;
}

template <class TAPosition>
void GlobalSubdivisionMultiGridRefiner<TAPosition>::nest_hierarchy()
{
	if(TAPosition::ValueType::Size == 1){
		UG_THROW("Error in GlobalSubdivisionMultiGridRefiner::project_to_subdivision_limit:\n"
				 "Currently only dimensions 2 and 3 are supported.\n");
	}

	ProjectHierarchyToSubdivisionLimit(*m_pMG, *m_pAPos);
}

template <class TAPosition>
void GlobalSubdivisionMultiGridRefiner<TAPosition>::refinement_step_ends()
{
	smooth();
}

template <class TAPosition>
void GlobalSubdivisionMultiGridRefiner<TAPosition>::smooth()
{
	if(TAPosition::ValueType::Size == 1){
		UG_THROW("Error in GlobalSubdivisionMultiGridRefiner::refinement_step_ends:\n"
				 "Currently only dimensions 2 and 3 are supported.\n");
	}

	ApplySmoothSubdivisionVolumesToTopLevel(*m_pMG, *m_pSH, *m_pMarkSH,
											*m_spLinearManifoldSH, m_bConstrained);
}


////////////////////////////////////////////////////////////////////////
//	Explicit instantiations
template class GlobalSubdivisionMultiGridRefiner<APosition1>;
template class GlobalSubdivisionMultiGridRefiner<APosition2>;
template class GlobalSubdivisionMultiGridRefiner<APosition>;


}//	end of namespace
