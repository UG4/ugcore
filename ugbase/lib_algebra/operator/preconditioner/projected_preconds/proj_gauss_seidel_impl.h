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

template<typename Matrix_type, typename Vector_type>
void forward_gs_step(Vector_type& c, const Matrix_type& A, const Vector_type& d,
		const size_t i, const number relaxFactor)
{
	typename Vector_type::value_type s;

	s = d[i];

	for(typename Matrix_type::const_row_iterator it = A.begin_row(i);
			it != A.end_row(i) && it.index() < i; ++it)
		// s -= it.value() * x[it.index()];
		MatMultAdd(s, 1.0, s, -1.0, it.value(), c[it.index()]);

	//	c[i] = relaxFactor * s / A(i,i)
	InverseMatMult(c[i], relaxFactor, A(i,i), s);
}

template<typename Matrix_type, typename Vector_type>
void backward_gs_step(Vector_type& c, const Matrix_type& A, const Vector_type& d,
		const size_t i, const number relaxFactor)
{
	typename Vector_type::value_type s;

	s = d[i];

	typename Matrix_type::const_row_iterator diag = A.get_connection(i, i);
	typename Matrix_type::const_row_iterator it = diag; ++it;

	for(; it != A.end_row(i); ++it)
		// s -= it.value() * x[it.index()];
		MatMultAdd(s, 1.0, s, -1.0, it.value(), c[it.index()]);

	// c[i] = relaxFactor * s/A(i,i)
	InverseMatMult(c[i], relaxFactor, diag.value(), s);
}


template <typename TAlgebra>
void
ProjGaussSeidel<TAlgebra>::
projected_precond_step(vector_type& c, vector_type& lastSol, const matrix_type& A, const vector_type& d)
{
	value_type tmpSol;

	for(size_t i = 0; i < c.size(); i++)
	{
		forward_gs_step(c, A, d, i, m_relax);

		//	compute temporary solution (solution of a common (forward) GaussSeidel-step)
		//	tmpSol := u_{s-1/2} = u_{s-1} + c
		tmpSol = lastSol[i] + c[i];

		//	perform a projection: check whether the temporary solution u_{s-1/2}
		//	fulfills the underlying constraint(s) or not
		if(m_spObsConstraint->lower_obs_set() && (!(m_spObsConstraint->upper_obs_set())))
		{
			//	only a lower obstacle is set
			m_spObsConstraint->correction_for_lower_obs(c, lastSol, i, tmpSol);
			continue;
		}

		if((!(m_spObsConstraint->lower_obs_set())) && m_spObsConstraint->upper_obs_set())
		{
			//	only an upper obstacle is set
			m_spObsConstraint->correction_for_upper_obs(c, lastSol, i, tmpSol);
			continue;
		}

		if(m_spObsConstraint->lower_obs_set() && m_spObsConstraint->upper_obs_set())
		{
			//	a lower and an upper obstacle are set
			m_spObsConstraint->correction_for_lower_and_upper_obs(c, lastSol, i, tmpSol);
			continue;
		}
	}
}

template <typename TAlgebra>
void
ProjBackwardGaussSeidel<TAlgebra>::
projected_precond_step(vector_type& c, vector_type& lastSol, const matrix_type& A, const vector_type& d)
{
	/*value_type tmpSol;

	if(c.size() == 0) return;
	size_t i = c.size()-1;
	do
	{
		backward_gs_step(c, A, d, i, m_relax);

		//	compute temporary solution (solution of a common backward GaussSeidel-step)
		//	tmpSol := u_{s-1/2} = u_{s-1} + c
		tmpSol = (*m_lastSol)[i] + c[i];

		//	perform a projection: check whether the temporary solution u_{s-1/2}
		//	fulfills the underlying constraint(s) or not
		if(m_bLowerObs && (!m_bUpperObs))
		{
			//	only a lower obstacle is set
			this->correction_for_lower_obs(c, i, tmpSol);
			continue;
		}

		if((!m_bLowerObs) && m_bUpperObs)
		{
			//	only an upper obstacle is set
			this->correction_for_upper_obs(c, i, tmpSol);
			continue;
		}

		if(m_bLowerObs && m_bUpperObs)
		{
			//	a lower and an upper obstacle are set
			this->correction_for_lower_and_upper_obs(c, i, tmpSol);
			continue;
		}

	} while(i-- != 0);*/

}

template <typename TAlgebra>
void
ProjSymmetricGaussSeidel<TAlgebra>::
projected_precond_step(vector_type& c, vector_type& lastSol, const matrix_type& A, const vector_type& d)
{
	/*value_type tmpSol;

	for(size_t i = 0; i < c.size(); i++)
	{
		//	1. perform a forward GaussSeidel step
		//	c1 = (D-L)^{-1} d
		forward_gs_step(c, A, d, i, m_relax);

		//	2. c2 = D c1
		MatMult(c[i], 1.0, A(i, i), c[i]);

		//	3. perform a backward GaussSeidel step
		//	c3 = (D-U)^{-1} c2
		backward_gs_step(c, A, c, i, m_relax);

		//	compute temporary solution (solution of a common symmetric GaussSeidel-step)
		//	tmpSol := u_{s-1/2} = u_{s-1} + c
		tmpSol = (*m_lastSol)[i] + c[i];

		//	perform a projection: check whether the temporary solution u_{s-1/2}
		//	fulfills the underlying constraint(s) or not
		if(m_bLowerObs && (!m_bUpperObs))
		{
			//	only a lower obstacle is set
			this->correction_for_lower_obs(c, i, tmpSol);
			continue;
		}

		if((!m_bLowerObs) && m_bUpperObs)
		{
			//	only an upper obstacle is set
			this->correction_for_upper_obs(c, i, tmpSol);
			continue;
		}

		if(m_bLowerObs && m_bUpperObs)
		{
			//	a lower and an upper obstacle are set
			this->correction_for_lower_and_upper_obs(c, i, tmpSol);
			continue;
		}
	}*/
}

} // end namespace ug

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__PROJ_GAUSS_SEIDEL_IMPL__ */
