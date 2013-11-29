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

template <typename TDomain, typename TAlgebra>
void
ScalarObstacle<TDomain, TAlgebra>::
correction_for_lower_obs(vector_type& c, vector_type& lastSol, const size_t index)
{
	//	TODO: in order to handle such cases, in which the obstacle is only set on some boundary-subset
	//	enter here something like
	//	if(obsSubset) do //oder if (indexIsInObsSubset)
	//	{	if (tmpSol - lowerObs) < 0.0 -> projection
	//		else lastSol = tmpSol;
	//	}
	//	else lastSol = tmpSol;

	//	compute unconstrained solution (solution of a common (forward) GaussSeidel-step)
	//	tmpSol := u_{s-1/2} = u_{s-1} + c
	value_type tmpSol = lastSol[index] + c[index];

	//	get index-th lower obstacle value
	const value_type& lowerObsVal = (*m_spVecOfLowObsValues)[index];

	if(GetSize(tmpSol) != GetSize(lowerObsVal))
		UG_THROW("size of tmpSol and size of upperObsVal need to be the same");

	for(size_t j = 0; j < GetSize(tmpSol); j++)
	{
		if ( (BlockRef(tmpSol, j) - BlockRef(lowerObsVal, j)) < 0.0)
		{
			//	u_{s-1/2} < lowerObsValue (:the lower constraint is not fulfilled)

			//	adjust correction c := u_s - u_{s-1} = m_obsVal - u_{s-1}
			BlockRef(c[index], j) = BlockRef(lowerObsVal, j) - BlockRef(lastSol[index], j);

			//	set new solution u_s to the obstacle value
			//	and store the current index in a vector for further treatment
			BlockRef(lastSol[index], j) = BlockRef(lowerObsVal, j);
			(*m_spLowerActiveInd).push_back(MultiIndex<2>(index, j) );
		}
		else
		{
			//	the 'tmpSol' is valid with respect to the lower constraints
			BlockRef(lastSol[index], j) = BlockRef(tmpSol, j);
		}
	}
}

template <typename TDomain, typename TAlgebra>
void
ScalarObstacle<TDomain, TAlgebra>::
correction_for_upper_obs(vector_type& c, vector_type& lastSol, const size_t index)
{
	//	compute unconstrained solution (solution of a common (forward) GaussSeidel-step)
	//	tmpSol := u_{s-1/2} = u_{s-1} + c
	value_type tmpSol = lastSol[index] + c[index];

	//	get index-th upper obstacle value
	const value_type& upperObsVal = (*m_spVecOfUpObsValues)[index];

	if(GetSize(tmpSol) != GetSize(upperObsVal))
		UG_THROW("size of tmpSol and size of upperObsVal need to be the same");

	for(size_t j = 0; j < GetSize(tmpSol); j++)
	{
		if ( (BlockRef(tmpSol, j) - BlockRef(upperObsVal, j)) > 0.0)
		{
			//	u_{s-1/2} > upperObsValue (:the upper constraint is not fulfilled)

			//	adjust correction c := u_s - u_{s-1} = m_obsVal - u_{s-1}
			BlockRef(c[index], j) = BlockRef(upperObsVal, j) - BlockRef(lastSol[index], j);

			//	set new solution u_s to the obstacle value
			//	and store the current index in a vector for further treatment
			BlockRef(lastSol[index], j) = BlockRef(upperObsVal, j);
			(*m_spUpperActiveInd).push_back(MultiIndex<2>(index, j) );
		}
		else
		{
			//	the 'tmpSol' is valid with respect to the upper constraints
			BlockRef(lastSol[index], j) = BlockRef(tmpSol, j);
		}
	}
}

template <typename TDomain, typename TAlgebra>
void
ScalarObstacle<TDomain, TAlgebra>::
correction_for_lower_and_upper_obs(vector_type& c, vector_type& lastSol, const size_t index)
{
	//	compute unconstrained solution (solution of a common (forward) GaussSeidel-step)
	//	tmpSol := u_{s-1/2} = u_{s-1} + c
	value_type tmpSol = lastSol[index] + c[index];

	//	get index-th lower obstacle value
	const value_type& upperObsVal = (*m_spVecOfUpObsValues)[index];
	const value_type& lowerObsVal = (*m_spVecOfLowObsValues)[index];

	if(GetSize(tmpSol) != GetSize(upperObsVal))
		UG_THROW("size of tmpSol and size of upperObsVal need to be the same");
	if(GetSize(tmpSol) != GetSize(lowerObsVal))
		UG_THROW("size of tmpSol and size of upperObsVal need to be the same");

	for(size_t j = 0; j < GetSize(tmpSol); j++)
	{
		if ( (BlockRef(tmpSol, j) - BlockRef(upperObsVal, j)) > 0.0)
		{
			//	u_{s-1/2} > upperObsValue (:the upper constraint is not fulfilled)

			//	adjust correction c := u_s - u_{s-1} = m_obsVal - u_{s-1}
			BlockRef(c[index], j)  = BlockRef(upperObsVal, j) - BlockRef(lastSol[index], j);

			//	set new solution u_s to the obstacle value
			//	and store the current index in a vector for further treatment
			BlockRef(lastSol[index], j) = BlockRef(upperObsVal, j);
			(*m_spUpperActiveInd).push_back(MultiIndex<2>(index, j) );
		}
		else
		{
			if ( (BlockRef(tmpSol, j) - BlockRef(lowerObsVal, j)) < 0.0)
			{
				//	u_{s-1/2} < lowerObsValue (:the lower constraint is not fulfilled)

				//	adjust correction c := u_s - u_{s-1} = m_obsVal - u_{s-1}
				BlockRef(c[index], j) = BlockRef(lowerObsVal, j) - BlockRef(lastSol[index], j);

				//	set new solution u_s to the obstacle value
				//	and store the current index in a vector for further treatment
				BlockRef(lastSol[index], j) = BlockRef(lowerObsVal, j);
				(*m_spLowerActiveInd).push_back(MultiIndex<2>(index, j) );
			}
			else
			{
				//	the 'tmpSol' is valid with respect to all constraints
				BlockRef(lastSol[index], j) = BlockRef(tmpSol, j);
			}
		}
	}
}

} // end namespace ug

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__SCALAR_OBSTACLE_IMPL__ */
