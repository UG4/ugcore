/*
 * assembled_linear_operator_impl.h
 *
 *  Created on: 24.02.2011
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR_IMPL__
#define __H__LIB_DISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR_IMPL__

#include "assembled_linear_operator.h"
#include "common/profiler/profiler.h"

#define PROFILE_ASS
#ifdef PROFILE_ASS
	#define ASS_PROFILE_FUNC()		PROFILE_FUNC()
	#define ASS_PROFILE_BEGIN(name)	PROFILE_BEGIN(name)
	#define ASS_PROFILE_END()		PROFILE_END()
#else
	#define ASS_PROFILE_FUNC()
	#define ASS_PROFILE_BEGIN(name)
	#define ASS_PROFILE_END()
#endif

namespace ug{

template <typename TDoFDistribution, typename TAlgebra>
bool
AssembledLinearOperator<TDoFDistribution, TAlgebra>::
init(const vector_type& u)
{
	if(m_pDoFDistribution == NULL)
	{
		UG_LOG("ERROR in AssembledLinearOperator::init: "
				"DoF Distribution not set.\n");
		return false;
	}

	if(m_pAss == NULL)
	{
		UG_LOG("ERROR in AssembledLinearOperator::init:"
				" Assembling routine not set.\n");
		return false;
	}

//	assemble matrix (depending on u, i.e. J(u))
	if(!m_pAss->assemble_jacobian(*this, u, *m_pDoFDistribution))
	{
		UG_LOG("ERROR in AssembledLinearOperator::init:"
				" Cannot assemble Jacobi matrix.\n");
		return false;
	}

//  we're done
	return true;
}

//	Initialize the operator
template <typename TDoFDistribution, typename TAlgebra>
bool
AssembledLinearOperator<TDoFDistribution, TAlgebra>::
init()
{
//	check if DoF Distribution is set
	if(m_pDoFDistribution == NULL)
	{
		UG_LOG("ERROR in AssembledLinearOperator::init:"
				" DoF Distribution not set.\n");
		return false;
	}
	if(m_pAss == NULL)
	{
		UG_LOG("ERROR in AssembledLinearOperator::init:"
				" Assembling routine not set.\n");
		return false;
	}

//	create vector dummy
	vector_type dummy; dummy.resize(m_pDoFDistribution->num_indices());

//	assemble only matrix
	if(!m_pAss->assemble_jacobian(*this, dummy, *m_pDoFDistribution))
	{
		UG_LOG("ERROR in AssembledLinearOperator::init:"
				" Cannot assemble Matrix.\n");
		return false;
	}

//	were done
	return true;
}

//	Initialize the operator
template <typename TDoFDistribution, typename TAlgebra>
bool
AssembledLinearOperator<TDoFDistribution, TAlgebra>::
init_op_and_rhs(vector_type& b)
{
//	todo: check that assembling is linear

//	check if DoF Distribution is set
	if(m_pDoFDistribution == NULL)
	{
		UG_LOG("ERROR in AssembledLinearOperator::init_op_and_rhs:"
				" DoF Distribution not set.\n");
		return false;
	}
	if(m_pAss == NULL)
	{
		UG_LOG("ERROR in AssembledLinearOperator::init_op_and_rhs:"
				" Assembling routine not set.\n");
		return false;
	}

//	resize rhs, since used as u dummy
	b.resize(m_pDoFDistribution->num_indices());

//	assemble matrix and rhs in one loop
	if(!m_pAss->assemble_linear(*this, b, b, *m_pDoFDistribution))
	{
		UG_LOG("ERROR in AssembledLinearOperator::init_op_and_rhs:"
				" Cannot assemble Matrix and Rhs.\n");
		return false;
	}

//	were done
	return true;
}

template <typename TDoFDistribution, typename TAlgebra>
bool
AssembledLinearOperator<TDoFDistribution, TAlgebra>::
apply(vector_type& d, const vector_type& c)
{
#ifdef UG_PARALLEL
	if(!c.has_storage_type(PST_CONSISTENT))
		{
			UG_LOG("WARNING: In 'AssembledLinearOperator::apply':"
					"Inadequate storage format of Vector c.");
			return false;
		}
#endif

//	perform check of sizes
	if(c.size() != this->num_cols() || d.size() != this->num_rows())
	{
		UG_LOG("ERROR in 'AssembledLinearOperator::apply': Size of matrix A ["<<
		        this->num_rows() << " x " << this->num_cols() << "] must match the "
		        "sizes of vectors x ["<<c.size()<<"], b ["<<d.size()<<"] for the "
		        " operation b = A*x. Maybe the operator is not initialized ?\n";)
		return false;
	}

//	Apply Matrix
	return base_type::apply(d, c);
}

//	Compute d := d - J(u)*c
template <typename TDoFDistribution, typename TAlgebra>
bool
AssembledLinearOperator<TDoFDistribution, TAlgebra>::
apply_sub(vector_type& d, const vector_type& c)
{
#ifdef UG_PARALLEL
	if(!d.has_storage_type(PST_ADDITIVE))
	{
		UG_LOG("ERROR: In 'AssembledLinearOperator::apply_sub':"
				"Inadequate storage format of Vector d.\n");
		return false;
	}
	if(!c.has_storage_type(PST_CONSISTENT))
	{
		UG_LOG("ERROR: In 'AssembledLinearOperator::apply_sub':"
				"Inadequate storage format of Vector c.\n");
		return false;
	}
#endif

//	check sizes
	if(c.size() != this->num_cols() || d.size() != this->num_rows())
	{
		UG_LOG("ERROR in AssembledLinearOperator::apply_sub: Size of matrix A ["<<
		        this->num_rows() << " x " << this->num_cols() << "] must match the "
		        "sizes of vectors x ["<<c.size()<<"], b ["<<d.size()<<"] for the "
		        " operation b -= A*x. Maybe the operator is not initialized ?\n";)
		return false;
	}

//	Apply Matrix
	return base_type::matmul_minus(d,c);
}


template <typename TDoFDistribution, typename TAlgebra>
bool
AssembledLinearOperator<TDoFDistribution, TAlgebra>::
set_dirichlet_values(vector_type& u)
{
//	checks
	if(m_pAss == NULL)
	{
		UG_LOG("ERROR in AssembledLinearOperator::set_dirichlet_values:"
				" Assembling routine not set.\n");
		return false;
	}
	if(m_pDoFDistribution == NULL)
	{
		UG_LOG("ERROR in AssembledLinearOperator::set_dirichlet_values:"
				" DoF Distribution not set.\n");
		return false;
	}

//	set dirichlet values etc.
	if(m_pAss->assemble_solution(u, *m_pDoFDistribution) != true)
	{
		UG_LOG("ERROR in AssembledLinearOperator::set_dirichlet_values:"
				" Cannot assemble solution.\n");
		return false;
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/// help function to assemble a linear operator
template <typename TDoFDistribution, typename TAlgebra>
bool AssembleLinearOperatorRhsAndSolution
		(AssembledLinearOperator<TDoFDistribution, TAlgebra>& op,
		 typename TAlgebra::vector_type& u,
		 typename TAlgebra::vector_type& b)
{
	ASS_PROFILE_BEGIN(ASS_AssembleLinearOperatorRhsAndSolution);

//	initialize operator
	ASS_PROFILE_BEGIN(ASS_InitOperatorAndRhs);
	if(!op.init_op_and_rhs(b))
	{
		UG_LOG("ERROR in 'AssembleLinearOperatorRhsAndSolution': Cannot init"
				" the operator (assembling failed).\n");
		return false;
	}
	ASS_PROFILE_END();

//	sets the dirichlet values in the solution
	ASS_PROFILE_BEGIN(ASS_SetDirValues);
	if(!op.set_dirichlet_values(u))
	{
		UG_LOG("ERROR in 'AssembleLinearOperatorRhsAndSolution': Cannot set"
				" the dirichlet values in the solution.\n");
		return false;
	}
	ASS_PROFILE_END();

	ASS_PROFILE_END();
//	done
	return true;
}


} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR_IMPL__ */
