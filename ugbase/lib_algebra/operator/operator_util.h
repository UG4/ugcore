/*
 * Copyright (c) 2010-2014:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef __H__LIB_ALGEBRA__OPERATOR__OPERATOR_UTIL__
#define __H__LIB_ALGEBRA__OPERATOR__OPERATOR_UTIL__


#include "common/profiler/profiler.h"
#include "common/util/smart_pointer.h"
#include "lib_algebra/operator/interface/linear_operator.h"
#include "lib_algebra/operator/interface/linear_operator_inverse.h"

namespace ug{

template <typename vector_type>
bool ApplyLinearSolver(	SmartPtr<ILinearOperator<vector_type> > A,
						vector_type& u, vector_type& b,
						SmartPtr<ILinearOperatorInverse<vector_type> > solver)
{
// step 1: Init Linear Inverse Operator
	PROFILE_BEGIN_GROUP(ALS_InitLinearSolver, "algebra");
	if(!solver->init(A))
	{
		UG_LOG("ApplyLinearSolver: Cannot init Inverse operator.\n");
		return false;
	}
	PROFILE_END_(ALS_InitLinearSolver);

// step 2: Apply Operator
	PROFILE_BEGIN(ALS_ApplyLinearSolver);
	if(!solver->apply_return_defect(u,b))
	{
		UG_LOG("ApplyLinearSolver: Cannot apply Inverse operator.\n");
		return false;
	}
	PROFILE_END_(ALS_ApplyLinearSolver);

//	done
	return true;
}


}
#endif
