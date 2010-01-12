/*
 * equation_impl.h
 *
 *  Created on: 17.11.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__EQUATION_IMPL__
#define __H__LIBDISCRETIZATION__EQUATION_IMPL__

namespace ug{

template <int d>
template <typename TElem>
bool Equation<d>::assemble_linear(TElem* elem, Matrix& mat, Vector& vec, NumericalSolution<d>& u)
{
	static DiscretizationScheme<TElem, d>* DiscScheme = &DiscretizationSchemes<TElem, d>::DiscretizationScheme(m_DiscretizationSchemeID);
	DiscScheme->reset_local_jacobian();
	DiscScheme->reset_local_defect();

	/* loop over all Differential Operators */
	for(uint i=0; i < m_ScalarDifferentialOperatorVector.size(); i++)
	{
		DiscScheme->add_op_to_local_jacobian(m_ScalarDifferentialOperatorVector[i], elem, u);
	}
	for(uint i=0; i < m_DivergenzDifferentialOperatorVector.size(); i++)
	{
		DiscScheme->add_op_to_local_jacobian(m_DivergenzDifferentialOperatorVector[i], elem, u);
	}
	if(DiscScheme->send_local_jacobian_to_global_jacobian(elem, mat, -1.0)==false) return false;

	/* loop over all RHS */
	for(uint i=0; i < m_RHSVector.size(); i++)
	{
		DiscScheme->add_rhs_to_local_defect(m_RHSVector[i], elem);
	}
	if(DiscScheme->send_local_defect_to_global_defect(elem, vec)==false) return false;

	return(true);
}

template <int d>
template <typename TElem>
bool Equation<d>::assemble_jacobian(TElem* elem, Matrix& mat, NumericalSolution<d>& u, number time, number s_m, number s_a)
{
	static DiscretizationScheme<TElem, d>* DiscScheme = &DiscretizationSchemes<TElem, d>::DiscretizationScheme(m_DiscretizationSchemeID);
	DiscScheme->reset_local_jacobian();

	/**** MASS MATRIX ****/
	/* loop over Time Operators */
	for(uint i=0; i < m_TimeOperatorVector.size(); i++)
	{
		DiscScheme->add_op_to_local_jacobian(m_TimeOperatorVector[i], elem, u);
	}
	/***** ADD LOCAL MATRIX TO GLOBAL MATRIX */
	if(DiscScheme->send_local_jacobian_to_global_jacobian(elem, mat, s_m)==false) return false;
	DiscScheme->reset_local_jacobian();

	/**** STIFFNESS MATRIX ****/
	/* loop over Scalar Operators */
	for(uint i=0; i < m_ScalarDifferentialOperatorVector.size(); i++)
	{
		DiscScheme->add_op_to_local_jacobian(m_ScalarDifferentialOperatorVector[i], elem, u);
	}

	/* loop over Differential Operators */
	for(uint i=0; i < m_DivergenzDifferentialOperatorVector.size(); i++)
	{
		DiscScheme->add_op_to_local_jacobian(m_DivergenzDifferentialOperatorVector[i], elem, u);
	}

	/***** ADD LOCAL MATRIX TO GLOBAL MATRIX */
	if(DiscScheme->send_local_jacobian_to_global_jacobian(elem, mat, s_a)==false) return false;

	return true;
}


template <int d>
template <typename TElem>
bool Equation<d>::assemble_defect(TElem* elem, Vector& vec, NumericalSolution<d>& u, number time, number s_m, number s_a)
{
	static DiscretizationScheme<TElem, d>* DiscScheme = &DiscretizationSchemes<TElem, d>::DiscretizationScheme(m_DiscretizationSchemeID);

	DiscScheme->reset_local_defect();

	/**** MASS MATRIX ****/
	/* loop over Time Operators */
	for(uint i=0; i < m_TimeOperatorVector.size(); i++)
	{
		DiscScheme->add_op_to_local_defect(m_TimeOperatorVector[i], elem, u);
	}
	/**** SEND LOCAL DEFECT TO GLOBAL DEFECT ****/
	if(DiscScheme->send_local_defect_to_global_defect(elem, vec, s_m)==false) return false;

	DiscScheme->reset_local_defect();

	/**** STIFFNESS MATRIX ****/
	/* loop over Scalar Operators */
	for(uint i=0; i < m_ScalarDifferentialOperatorVector.size(); i++)
	{
		DiscScheme->add_op_to_local_defect(m_ScalarDifferentialOperatorVector[i], elem, u);
	}

	/* loop over Differential Operators */
	for(uint i=0; i < m_DivergenzDifferentialOperatorVector.size(); i++)
	{
		DiscScheme->add_op_to_local_defect(m_DivergenzDifferentialOperatorVector[i], elem, u);
	}

	/**** RHS ****/
	/* loop over all RHS */
	for(uint i=0; i < m_RHSVector.size(); i++)
	{
		DiscScheme->add_rhs_to_local_defect(m_RHSVector[i], elem);
	}

	/**** SEND LOCAL DEFECT TO GLOBAL DEFECT ****/
	if(DiscScheme->send_local_defect_to_global_defect(elem, vec, s_a)==false) return false;

	return true;
}



template <int d>
template <typename TElem>
bool Equation<d>::assemble_defect(TElem* elem, Vector& vec, NumericalSolution<d>& u)
{
	static DiscretizationScheme<TElem, d>* DiscScheme = &DiscretizationSchemes<TElem, d>::DiscretizationScheme(m_DiscretizationSchemeID);

	DiscScheme->reset_local_defect();

	/**** STIFFNESS MATRIX ****/
	/* loop over Scalar Operators */
	for(uint i=0; i < m_ScalarDifferentialOperatorVector.size(); i++)
	{
		DiscScheme->add_op_to_local_defect(m_ScalarDifferentialOperatorVector[i], elem, u);
	}

	/* loop over Differential Operators */
	for(uint i=0; i < m_DivergenzDifferentialOperatorVector.size(); i++)
	{
		DiscScheme->add_op_to_local_defect(m_DivergenzDifferentialOperatorVector[i], elem, u);
	}

	/**** RHS ****/
	/* loop over all RHS */
	for(uint i=0; i < m_RHSVector.size(); i++)
	{
		DiscScheme->add_rhs_to_local_defect(m_RHSVector[i], elem);
	}

	/**** SEND LOCAL DEFECT TO GLOBAL DEFECT ****/
	if(DiscScheme->send_local_defect_to_global_defect(elem, vec)==false) return false;

	return true;
}


template <int d>
template <typename TElem>
bool Equation<d>::assemble_jacobian(TElem* elem, Matrix& mat, NumericalSolution<d>& u)
{
	static DiscretizationScheme<TElem, d>* DiscScheme = &DiscretizationSchemes<TElem, d>::DiscretizationScheme(m_DiscretizationSchemeID);

	DiscScheme->reset_local_jacobian();

	/**** STIFFNESS MATRIX ****/
	/* loop over Scalar Operators */
	for(uint i=0; i < m_ScalarDifferentialOperatorVector.size(); i++)
	{
		DiscScheme->add_op_to_local_jacobian(m_ScalarDifferentialOperatorVector[i], elem, u);
	}

	/* loop over Differential Operators */
	for(uint i=0; i < m_DivergenzDifferentialOperatorVector.size(); i++)
	{
		DiscScheme->add_op_to_local_jacobian(m_DivergenzDifferentialOperatorVector[i], elem, u);
	}

	/***** ADD LOCAL MATRIX TO GLOBAL MATRIX */
	if(DiscScheme->send_local_jacobian_to_global_jacobian(elem, mat)==false) return false;

	return true;
}


}
#endif /* __H__LIBDISCRETIZATION__EQUATION_IMPL__ */
