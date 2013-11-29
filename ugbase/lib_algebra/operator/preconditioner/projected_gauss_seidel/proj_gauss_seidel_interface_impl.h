/*
 * proj_gauss_seidel_interface_impl.h
 *
 *  Created on: 13.11.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL_INTERFACE_IMPL__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL_INTERFACE_IMPL__

#include "proj_gauss_seidel_interface.h"

#define PROFILE_PROJ_GAUSS_SEIDEL
#ifdef PROFILE_PROJ_GAUSS_SEIDEL
	#define PROJ_GAUSS_SEIDEL_PROFILE_FUNC()		PROFILE_FUNC()
	#define PROJ_GAUSS_SEIDEL_PROFILE_BEGIN(name)	PROFILE_BEGIN_GROUP(name, "IProjGaussSeidel")
	#define PROJ_GAUSS_SEIDEL_PROFILE_END()		PROFILE_END()
#else
	#define PROJ_GAUSS_SEIDEL_PROFILE_FUNC()
	#define PROJ_GAUSS_SEIDEL_PROFILE_BEGIN(name)
	#define PROJ_GAUSS_SEIDEL_PROFILE_END()
#endif

namespace ug{

template <typename TAlgebra>
bool
IProjGaussSeidel<TAlgebra>::
init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u)
{
	PROFILE_FUNC_GROUP("IProjGaussSeidel");

	base_type::init(J,u);

//	remember solution
	m_spSol = u.clone();

//	remember, that operator has been initialized
	m_bInit = true;

//	init the obstacle constraint class
	if(m_bObsCons)
		m_spObsConstraint->init(u);

//	(ugly) hint, that usual damping (x += damp * c) does not make sense for the projected
//	GaussSeidel-method.
/*	const number kappa = this->damping()->damping(u, u, m_spOperator);
	if(kappa != 1.0){
		UG_THROW("IProjGaussSeidel::set_damp': Ususal damping is not possible "
				"for IProjGaussSeidel! Use 'set_sor_relax' instead!");
	}*/

	return true;
}

template <typename TAlgebra>
void
IProjGaussSeidel<TAlgebra>::project_correction(vector_type& c, const size_t i)
{
	//	compute unconstrained solution (solution of a common (forward) GaussSeidel-step)
	//	tmpSol := u_{s-1/2} = u_{s-1} + c
	//value_type tmpSol = (*m_spSol)[i] + c[i];

	//TODO: check, if index i is in m_vAlgIndicesOfObsSubsets
	//if (m_spObsConstraint->nrOfIndicesOfObsSubsets() != 0.0)
	//if (m_spObsConstraint->m_vAlgIndicesOfObsSubsets.size() != 0.0)
	//	if (!m_spObsConstraint->index_is_in_obs_subset(i))
	//		return;

	//	perform a projection: check whether the temporary solution u_{s-1/2}
	//	fulfills the underlying constraint(s) or not
	if(m_spObsConstraint->lower_obs_set() && (!(m_spObsConstraint->upper_obs_set())))
	{
		//	only a lower obstacle is set
		m_spObsConstraint->correction_for_lower_obs(c, (*m_spSol), i);
	}

	if((!(m_spObsConstraint->lower_obs_set())) && m_spObsConstraint->upper_obs_set())
	{
		//	only an upper obstacle is set
		m_spObsConstraint->correction_for_upper_obs(c, (*m_spSol), i);
	}

	if(m_spObsConstraint->lower_obs_set() && m_spObsConstraint->upper_obs_set())
	{
		//	a lower and an upper obstacle are set
		m_spObsConstraint->correction_for_lower_and_upper_obs(c, (*m_spSol), i);
	}
}

template <typename TAlgebra>
bool
IProjGaussSeidel<TAlgebra>::
apply(vector_type &c, const vector_type& d)
{
	PROFILE_FUNC_GROUP("IProjGaussSeidel");
//	Check that operator is initialized
	if(!m_bInit)
	{
		UG_LOG("ERROR in 'IProjGaussSeidel::apply': Iterator not initialized.\n");
		return false;
	}

	//	reset active indices
	m_spObsConstraint->reset_active_indices();

	base_type::apply(c, d);

	return true;
}

template <typename TAlgebra>
void
IProjGaussSeidel<TAlgebra>::
adjust_defect(vector_type& d)
{
	for (std::vector<MultiIndex<2> >::iterator itActiveInd = (*(m_spObsConstraint->lower_active_indices())).begin();
					itActiveInd < (*(m_spObsConstraint->lower_active_indices())).end(); ++itActiveInd)
	{
		//	check, if Ax <= b. For that case the new defect is set to zero,
		//	since all equations/constraints are fulfilled
		if (BlockRef(d[(*itActiveInd)[0]], (*itActiveInd)[1]) < 0.0)
			BlockRef(d[(*itActiveInd)[0]], (*itActiveInd)[1]) = 0.0;
	}

	for (std::vector<MultiIndex<2> >::iterator itActiveInd = (*(m_spObsConstraint->upper_active_indices())).begin();
				itActiveInd < (*(m_spObsConstraint->upper_active_indices())).end(); ++itActiveInd)
	{
		//	check, if Ax >= b. For that case the new defect is set to zero,
		//	since all equations/constraints are fulfilled
		if (BlockRef(d[(*itActiveInd)[0]], (*itActiveInd)[1]) > 0.0)
			BlockRef(d[(*itActiveInd)[0]], (*itActiveInd)[1]) = 0.0;
	}
}

template <typename TAlgebra>
bool
IProjGaussSeidel<TAlgebra>::
apply_update_defect(vector_type &c, vector_type& d)
{
	PROFILE_FUNC_GROUP("IProjGaussSeidel");

	base_type::apply_update_defect(c, d);

	//	adjust defect of the active indices for the case that a constraint is set
	if(m_bObsCons)
		adjust_defect(d);

	return true;
}


} // end namespace ug

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL_INTERFACE_IMPL__ */
