/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR_IMPL__
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR_IMPL__

#include "assembled_non_linear_operator.h"

namespace ug{

//	Prepare functions
template <typename TAlgebra>
void
AssembledOperator<TAlgebra>::prepare(vector_type& uIn)
{
	PROFILE_BEGIN_GROUP(AssembledOperator_prepare, "discretization");
	if(m_spAss.invalid())
		UG_THROW("Discretization not set.");

// 	Set Dirichlet - Nodes to exact values (any constraints in general)
	try{
		m_spAss->adjust_solution(uIn, m_gridLevel);
	}
	UG_CATCH_THROW("Cannot set constraints in solution.");
}

// 	Compute d = L(u)
template <typename TAlgebra>
void
AssembledOperator<TAlgebra>::apply(vector_type& dOut, const vector_type& uIn)
{
	PROFILE_BEGIN_GROUP(AssembledOperator_apply, "discretization");
	if(m_spAss.invalid())
		UG_THROW("Discretization not set.");

//  assemble defect
	try{
		m_spAss->assemble_defect(dOut, uIn, m_gridLevel);
	}
	UG_CATCH_THROW("Could not assemble defect. Aborting.");
}

} // end namepace ug

#endif /*__H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR_IMPL__*/
