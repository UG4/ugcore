/*
 * Copyright (c) 2020:  G-CSC, Goethe University Frankfurt
 * Author: Lukas Larisch
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

#ifndef __UG__LIB_DISC__ORDERING_STRATEGIES_EXECUTION_DOFS_ORDERING__
#define __UG__LIB_DISC__ORDERING_STRATEGIES_EXECUTION_DOFS_ORDERING__

#include <vector>

#include "lib_disc/function_spaces/approximation_space.h"

#include "../../../lib_algebra/ordering_strategies/execution/IExecuteOrdering.h"
#include "../../../lib_algebra/ordering_strategies/execution/util.cpp"

#include "../../../common/code_marker.h" //error()

namespace ug{

/*
	Orders all dof_distribution(index) according to an ordering. dof_distribution triggers reordering
	of handled matrices.
*/

//TODO use GridLevel
template <typename TDomain, typename O_t=std::vector<size_t> >
class DOFsOrdering : public IExecuteOrdering{
public:
	DOFsOrdering(SmartPtr<ug::ApproximationSpace<TDomain> > approx_in, O_t &o_in, unsigned idx_in)
		: m_spApprox(approx_in), o(o_in), idx(idx_in){}
	virtual ~DOFsOrdering(){}

	void execute(){
		if(!is_permutation(o)){
			std::cerr << "[DOFsOrdering::execute] Not a permutation!" << std::endl;
			error();
		}

		m_spApprox->dof_distributions()[idx]->permute_indices(o);
	}

private:
	SmartPtr<ug::ApproximationSpace<TDomain> > m_spApprox;
	O_t &o;
	size_t idx;
};

} //namespace

#endif //guard
