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
init()
{
	base_type::initObsValues();

	base_type::get_dof_on_obs_subset(m_vLowerObsDoFs);
	base_type::get_obstacle_value_map(m_mScalarLowerObsValues);
}

template <typename TDomain, typename TAlgebra>
void
ScalarLowerObstacle<TDomain, TAlgebra>::
project_on_admissible_set(value_type& c_i, value_type& sol_i, const size_t i)
{
	//	loop all components of the correction at index i 'c_i'
	for(size_t comp = 0; comp < GetSize(c_i); comp++)
	{
		DoFIndex dof = DoFIndex(i, comp);

		//	check, if the dof lies in an obstacle subset: if not -> continue!
		if (!base_type::dof_lies_on_obs_subset(dof, m_vLowerObsDoFs))
			continue;

		//	tmpSol := u_{s-1/2} = u_{s-1} + c
		const number tmpSol = BlockRef(sol_i, comp) + BlockRef(c_i, comp);

		//	get lower obstacle value corresponding to the dof
		const number lowerObsVal = m_mScalarLowerObsValues[dof];

		//	check, if dof is admissible
		if ((tmpSol - lowerObsVal) < 0.0)
		{
			//	not admissible -> active DoF
			m_vActiveDofs.push_back(dof);

			BlockRef(c_i, comp) = lowerObsVal - BlockRef(sol_i, comp);
			BlockRef(sol_i, comp) = lowerObsVal;
		}
		else{
			// dof is admissible -> do regular update
			BlockRef(sol_i, comp) += BlockRef(c_i, comp);
		}
	} //end(for(comp))
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
		if (BlockRef(d[(*itActiveInd)[0]], (*itActiveInd)[1]) < 0.0)
			BlockRef(d[(*itActiveInd)[0]], (*itActiveInd)[1]) = 0.0;
	}
}


//////////////////////////////
//	SCALAR UPPER OBSTACLE
//////////////////////////////

template <typename TDomain, typename TAlgebra>
void
ScalarUpperObstacle<TDomain, TAlgebra>::
init()
{
	base_type::initObsValues();

	base_type::get_dof_on_obs_subset(m_vUpperObsDoFs);
	base_type::get_obstacle_value_map(m_mScalarUpperObsValues);
}


template <typename TDomain, typename TAlgebra>
void
ScalarUpperObstacle<TDomain, TAlgebra>::
project_on_admissible_set(value_type& c_i, value_type& sol_i, const size_t i)
{
	//	loop all components of the correction at index i 'c_i'
	for(size_t comp = 0; comp < GetSize(c_i); comp++)
	{
		DoFIndex dof = DoFIndex(i, comp);

		//	check, if the dof lies in an obstacle subset: if not -> continue!
		if (!base_type::dof_lies_on_obs_subset(dof, m_vUpperObsDoFs))
			continue;

		//	tmpSol := u_{s-1/2} = u_{s-1} + c
		const number tmpSol = BlockRef(sol_i, comp) + BlockRef(c_i, comp);

		//	get upper obstacle value corresponding to the dof
		const number upperObsVal = m_mScalarUpperObsValues[dof];

		//	check, if dof is admissible
		if ((tmpSol - upperObsVal) > 0.0)
		{
			//	not admissible -> active DoF
			m_vActiveDofs.push_back(dof);

			BlockRef(c_i, comp) = upperObsVal - BlockRef(sol_i, comp);
			BlockRef(sol_i, comp) = upperObsVal;
		}
		else{
			// dof is admissible -> do regular update
			BlockRef(sol_i, comp) += BlockRef(c_i, comp);
		}
	} //end(for(comp))
}

template <typename TDomain, typename TAlgebra>
void
ScalarUpperObstacle<TDomain, TAlgebra>::
adjust_defect(vector_type& d)
{
	for (std::vector<MultiIndex<2> >::iterator itActiveInd = m_vActiveDofs.begin();
			itActiveInd < m_vActiveDofs.end(); ++itActiveInd)
	{
		//	check, if Ax <= b. For that case the new defect is set to zero,
		//	since all equations/constraints are fulfilled
		if (BlockRef(d[(*itActiveInd)[0]], (*itActiveInd)[1]) > 0.0)
			BlockRef(d[(*itActiveInd)[0]], (*itActiveInd)[1]) = 0.0;
	}
}

} // end namespace ug

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__SCALAR_OBSTACLE_IMPL__ */
