/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Raphael Prohl
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
adjust_defect_to_constraint(vector_type& d)
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
adjust_defect_to_constraint(vector_type& d)
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
