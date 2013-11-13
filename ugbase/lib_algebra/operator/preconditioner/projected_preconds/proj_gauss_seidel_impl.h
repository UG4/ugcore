/*
 * proj_gauss_seidel_impl.h
 *
 *  Created on: 10.10.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__PROJ_GAUSS_SEIDEL_IMPL__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__PROJ_GAUSS_SEIDEL_IMPL__

#include "proj_gauss_seidel.h"

namespace ug{

template <typename TAlgebra>
void
ProjGaussSeidel<TAlgebra>::
projected_precond_step(vector_type& c, const matrix_type& A, const vector_type& d)
{
	//	TODO: matrix_type correct or MatrixOperator<matrix_type, vector_type>?
	value_type s, tmpSol;
	const number relaxFactor = m_relax;

	for(size_t i = 0; i < c.size(); i++)
	{
		s = d[i];

		for(typename matrix_type::const_row_iterator it = A.begin_row(i);
				it != A.end_row(i) && it.index() < i; ++it)
			// s -= it.value() * x[it.index()];
			MatMultAdd(s, 1.0, s, -1.0, it.value(), c[it.index()]);

		//	c[i] = relaxFactor * s / A(i,i)
		InverseMatMult(c[i], relaxFactor, A(i,i), s);

		//	compute temporary solution (solution of a common (forward) GaussSeidel-step)
		//	tmpSol := u_{s-1/2} = u_{s-1} + c
		tmpSol = m_lastSol[i] + c[i];

		//	perform a projection: check whether the temporary solution u_{s-1/2}
		//	fulfills the underlying constraint(s) or not
		if(m_bLowerObs && (!m_bUpperObs))
		{
			//	only a lower obstacle is set
			correction_for_lower_obs(c, i, tmpSol);
			continue;
		}

		if((!m_bLowerObs) && m_bUpperObs)
		{
			//	only an upper obstacle is set
			correction_for_upper_obs(c, i, tmpSol);
			continue;
		}

		if(m_bLowerObs && m_bUpperObs)
		{
			//	a lower and an upper obstacle are set
			correction_for_lower_and_upper_obs(c, i, tmpSol);
			continue;
		}

		//UG_LOG("m_lastSol[" << i << "]: " << m_lastSol[i] << "\n");
	}
}

} // end namespace ug

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__PROJ_GAUSS_SEIDEL_IMPL__ */
