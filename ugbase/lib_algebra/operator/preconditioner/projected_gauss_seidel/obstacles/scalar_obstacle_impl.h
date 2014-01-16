/*
 * scalar_obstacle_impl.h
 *
 *  Created on: 25.11.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__SCALAR_OBSTACLE_IMPL__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__SCALAR_OBSTACLE_IMPL__

#include "scalar_obstacle.h"

namespace ug{

//////////////////////////////
//	SCALAR LOWER OBSTACLE
//////////////////////////////

template <typename TDomain, typename TAlgebra>
void
ScalarLowerObstacle<TDomain, TAlgebra>::
adjust_sol_and_cor(value_type& sol_i, value_type& c_i, bool& dofIsActive,
		const DoFIndex& dof)
{
	const size_t comp = dof[1];

	//	tmpSol := u_{s-1/2} = u_{s-1} + c
	const number tmpSol = BlockRef(sol_i, comp) + BlockRef(c_i, comp);

	//	get lower obstacle value corresponding to the dof
	const number obsVal = m_mObstacleValues[dof];

	//	check, if dof is active (tmpSol <= obsVal)
	if (!(tmpSol > obsVal))
	{
		//	is active DoF
		m_vActiveDofs.push_back(dof);

		//	adjust correction & set solution to obstacle-value
		BlockRef(c_i, comp) = obsVal - BlockRef(sol_i, comp);
		BlockRef(sol_i, comp) = obsVal;
		dofIsActive = true;
	}
}

template <typename TDomain, typename TAlgebra>
void
ScalarLowerObstacle<TDomain, TAlgebra>::
adjust_defect(vector_type& d)
{
	for (std::vector<MultiIndex<2> >::iterator itActiveInd = m_vActiveDofs.begin();
			itActiveInd < m_vActiveDofs.end(); ++itActiveInd)
	{
		//	check, if Ax <= b. For that case the new defect is set to zero,
		//	since all equations/constraints are fulfilled
		number defect = BlockRef(d[(*itActiveInd)[0]], (*itActiveInd)[1]);
		if (defect < 0.0)
			BlockRef(d[(*itActiveInd)[0]], (*itActiveInd)[1]) = 0.0;
	}
}

template <typename TDomain, typename TAlgebra>
void
ScalarLowerObstacle<TDomain, TAlgebra>::
restrict_obs_values()
{}

//////////////////////////////
//	SCALAR UPPER OBSTACLE
//////////////////////////////

template <typename TDomain, typename TAlgebra>
void
ScalarUpperObstacle<TDomain, TAlgebra>::
adjust_sol_and_cor(value_type& sol_i, value_type& c_i, bool& dofIsActive,
		const DoFIndex& dof)
{
	const size_t comp = dof[1];

	//	tmpSol := u_{s-1/2} = u_{s-1} + c
	const number tmpSol = BlockRef(sol_i, comp) + BlockRef(c_i, comp);

	//	get upper obstacle value corresponding to the dof
	const number obsVal = m_mObstacleValues[dof];

	//	check, if dof is active (tmpSol >= obsVal)
	if (!(tmpSol < obsVal))
	{
		//	is active DoF
		m_vActiveDofs.push_back(dof);

		//	adjust correction & set solution to obstacle-value
		BlockRef(c_i, comp) = obsVal - BlockRef(sol_i, comp);
		BlockRef(sol_i, comp) = obsVal;
		dofIsActive = true;
	}
	//UG_LOG("dof " <<dof<< " is active: " <<dofIsActive<<"\n");
}

template <typename TDomain, typename TAlgebra>
void
ScalarUpperObstacle<TDomain, TAlgebra>::
adjust_defect(vector_type& d)
{
	for (std::vector<MultiIndex<2> >::iterator itActiveInd = m_vActiveDofs.begin();
			itActiveInd < m_vActiveDofs.end(); ++itActiveInd)
	{
		//	check, if Ax > b. For that case the new defect is set to zero,
		//	since all equations/constraints are fulfilled
		//UG_LOG("adjust_defect: " << (*itActiveInd)[0] <<","<< (*itActiveInd)[1] << "\n");
		number defect = BlockRef(d[(*itActiveInd)[0]], (*itActiveInd)[1]);
		if (defect > 0.0)
		{
			//UG_LOG("defect > 0 \n");
			BlockRef(d[(*itActiveInd)[0]], (*itActiveInd)[1]) = 0.0;
		}
	}
}

template <typename TDomain, typename TAlgebra>
void
ScalarUpperObstacle<TDomain, TAlgebra>::
restrict_obs_values()
{}

} // end namespace ug

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__SCALAR_OBSTACLE_IMPL__ */
