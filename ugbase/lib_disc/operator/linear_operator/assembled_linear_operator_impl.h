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

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR_IMPL__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR_IMPL__

#include "assembled_linear_operator.h"
#include "common/profiler/profiler.h"

#define PROFILE_ASS
#ifdef PROFILE_ASS
	#define ASS_PROFILE_FUNC()		PROFILE_FUNC()
	#define ASS_PROFILE_BEGIN(name)	PROFILE_BEGIN_GROUP(name, "discretization")
	#define ASS_PROFILE_END()		PROFILE_END()
#else
	#define ASS_PROFILE_FUNC()
	#define ASS_PROFILE_BEGIN(name)
	#define ASS_PROFILE_END()
#endif

namespace ug{

template <typename TAlgebra>
void
AssembledLinearOperator<TAlgebra>::init(const vector_type& u)
{
	if(m_spAss.invalid())
		UG_THROW("AssembledLinearOperator: Assembling routine not set.");

//	assemble matrix (depending on u, i.e. J(u))
	try{
		m_spAss->assemble_jacobian(*this, u, m_gridLevel);
	}
	UG_CATCH_THROW("AssembledLinearOperator: Cannot assemble Jacobi matrix.");
}

//	Initialize the operator
template <typename TAlgebra>
void
AssembledLinearOperator<TAlgebra>::init()
{
	if(m_spAss.invalid())
		UG_THROW("AssembledLinearOperator: Assembling routine not set.");

//	create vector dummy
	vector_type dummy;

//	assemble only matrix
	try{
		m_spAss->assemble_linear(*this, dummy, m_gridLevel);
	}
	UG_CATCH_THROW("AssembledLinearOperator::init: Cannot assemble Matrix.");
}

//	Initialize the operator
template <typename TAlgebra>
void
AssembledLinearOperator<TAlgebra>::init_op_and_rhs(vector_type& b)
{
//	todo: check that assembling is linear

	if(m_spAss.invalid())
		UG_THROW("AssembledLinearOperator: Assembling routine not set.");

//	assemble matrix and rhs in one loop
	try{
		m_spAss->assemble_linear(*this, b, m_gridLevel);
	}
	UG_CATCH_THROW("AssembledLinearOperator::init_op_and_rhs:"
						" Cannot assemble Matrix and Rhs.");
}

template <typename TAlgebra>
void
AssembledLinearOperator<TAlgebra>::apply(vector_type& d, const vector_type& c)
{
#ifdef UG_PARALLEL
	if(!c.has_storage_type(PST_CONSISTENT))
		UG_THROW("Inadequate storage format of Vector c.");
#endif

//	perform check of sizes
	if(c.size() != this->num_cols() || d.size() != this->num_rows())
		UG_THROW("AssembledLinearOperator::apply: Size of matrix A ["<<
		        this->num_rows() << " x " << this->num_cols() << "] must match the "
		        "sizes of vectors x ["<<c.size()<<"], b ["<<d.size()<<"] for the "
		        " operation b = A*x. Maybe the operator is not initialized ?");

//	Apply Matrix
	base_type::apply(d, c);
}

//	Compute d := d - J(u)*c
template <typename TAlgebra>
void
AssembledLinearOperator<TAlgebra>::apply_sub(vector_type& d, const vector_type& c)
{
#ifdef UG_PARALLEL
	if(!d.has_storage_type(PST_ADDITIVE))
		UG_THROW("Inadequate storage format of Vector d.");
	if(!c.has_storage_type(PST_CONSISTENT))
		UG_THROW("Inadequate storage format of Vector c.");
#endif

//	check sizes
	if(c.size() != this->num_cols() || d.size() != this->num_rows())
		UG_THROW("AssembledLinearOperator::apply_sub: Size of matrix A ["<<
		        this->num_rows() << " x " << this->num_cols() << "] must match the "
		        "sizes of vectors x ["<<c.size()<<"], b ["<<d.size()<<"] for the "
		        " operation b -= A*x. Maybe the operator is not initialized ?");

//	Apply Matrix
	//UG_LOG_ALL_PROCS("assembled linear iterator apply sub")
	base_type::matmul_minus(d,c);
}


template <typename TAlgebra>
void AssembledLinearOperator<TAlgebra>::set_dirichlet_values(vector_type& u)
{
//	checks
	if(m_spAss.invalid())
		UG_THROW("AssembledLinearOperator: Assembling routine not set.");

//	set dirichlet values etc.
	try{
		m_spAss->adjust_solution(u, m_gridLevel);
	}
	UG_CATCH_THROW("AssembledLinearOperator::set_dirichlet_values:"
				" Cannot assemble solution.");
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/// help function to assemble a linear operator
template <typename TAlgebra>
void AssembleLinearOperatorRhsAndSolution
		(AssembledLinearOperator<TAlgebra>& op,
		 typename TAlgebra::vector_type& u,
		 typename TAlgebra::vector_type& b)
{
	ASS_PROFILE_BEGIN(ASS_AssembleLinearOperatorRhsAndSolution);

//	initialize operator
	ASS_PROFILE_BEGIN(ASS_InitOperatorAndRhs);
	try{
		op.init_op_and_rhs(b);
	}UG_CATCH_THROW("Cannot init the operator (assembling failed).");
	ASS_PROFILE_END();

//	sets the dirichlet values in the solution
	ASS_PROFILE_BEGIN(ASS_SetDirValues);
	try{
		op.set_dirichlet_values(u);
	} UG_CATCH_THROW("Cannot set the dirichlet values in the solution.");
	ASS_PROFILE_END();

	ASS_PROFILE_END();
}


} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR_IMPL__ */
