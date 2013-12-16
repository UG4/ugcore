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

template <typename TDomain, typename TAlgebra>
bool
IProjGaussSeidel<TDomain,TAlgebra>::
init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u)
{
	PROFILE_FUNC_GROUP("IProjGaussSeidel");

//	call GaussSeidelBase-init
	base_type::init(J,u);

//	remember solution
	m_spSol = u.clone();

//	remember, that operator has been initialized
	m_bInit = true;

//	(ugly) hint, that usual damping (x += damp * c) does not make sense for the projected
//	GaussSeidel-method.
/*	const number kappa = this->damping()->damping(u, u, m_spOperator);
	if(kappa != 1.0){
		UG_THROW("IProjGaussSeidel::set_damp': Ususal damping is not possible "
				"for IProjGaussSeidel! Use 'set_sor_relax' instead!");
	}*/

	return true;
}

template <typename TDomain, typename TAlgebra>
void
IProjGaussSeidel<TDomain,TAlgebra>::
project_correction(value_type& c_i, const size_t i)
{
	if(!m_bObsCons)
		return;

	typedef typename vector<SmartPtr<IObstacleConstraint<TDomain,TAlgebra> > >::iterator iter_type;
	iter_type iterEnd = m_spvObsConstraint.end();

	for(size_t comp = 0; comp < GetSize(c_i); comp++)
	{
		DoFIndex dof = DoFIndex(i, comp);

		//	loop all obstacle constraint, which are set
		//	& perform a projection: check whether the temporary solution u_{s-1/2}
		//	fulfills the underlying constraint(s) or not
		bool dofIsAdmissible = true;
		bool dofIsObsDoF = false;

		//	set iterator to the first obstacle constraint
		iter_type iter = m_spvObsConstraint.begin();
		for( ; iter != iterEnd; iter++)
		{
			//	check, if the dof lies in an obstacle subset: if not -> continue!
			if (!((*iter)->dof_lies_in_obs_subset(dof)))
				continue;

			dofIsObsDoF = true;

			(*iter)->adjust_sol_and_cor((*m_spSol)[i], c_i, dofIsAdmissible, dof);
		}

		if (dofIsObsDoF && dofIsAdmissible)
		{
			// dof is admissible -> do regular solution update
			BlockRef((*m_spSol)[i], comp) += BlockRef(c_i, comp);
		}
	}
}

template <typename TDomain, typename TAlgebra>
bool
IProjGaussSeidel<TDomain,TAlgebra>::
apply(vector_type &c, const vector_type& d)
{
	PROFILE_FUNC_GROUP("IProjGaussSeidel");
//	Check that operator is initialized
	if(!m_bInit)
	{
		UG_LOG("ERROR in 'IProjGaussSeidel::apply': Iterator not initialized.\n");
		return false;
	}

	//	loop all obstacle constraints, which are set & reset its active dofs
	if(m_bObsCons)
	{
		typedef typename vector<SmartPtr<IObstacleConstraint<TDomain,TAlgebra> > >::iterator iter_type;
		iter_type iter = m_spvObsConstraint.begin();
		iter_type iterEnd = m_spvObsConstraint.end();

		for( ; iter != iterEnd; iter++)
			(*iter)->reset_active_dofs();
	}

	base_type::apply(c, d);

	return true;
}

template <typename TDomain, typename TAlgebra>
bool
IProjGaussSeidel<TDomain,TAlgebra>::
apply_update_defect(vector_type &c, vector_type& d)
{
	PROFILE_FUNC_GROUP("IProjGaussSeidel");

	base_type::apply_update_defect(c, d);

	//	adjust defect for the active dofs
	if(m_bObsCons)
	{
		typedef typename vector<SmartPtr<IObstacleConstraint<TDomain,TAlgebra> > >::iterator iter_type;
		iter_type iter = m_spvObsConstraint.begin();
		iter_type iterEnd = m_spvObsConstraint.end();

		for( ; iter != iterEnd; iter++)
			(*iter)->adjust_defect(d);
	}

	return true;
}


} // end namespace ug

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL_INTERFACE_IMPL__ */
