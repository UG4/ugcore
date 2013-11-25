/*
 * scalar_obstacle_impl.h
 *
 *  Created on: 25.11.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__SCALAR_OBSTACLE_IMPL__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__SCALAR_OBSTACLE_IMPL__

#include "scalar_obstacle.h"

namespace ug{

template <typename TAlgebra>
void
ScalarObstacle<TAlgebra>::
correction_for_lower_obs(vector_type& c, vector_type& lastSol, const size_t index, const value_type& tmpSol)
{
	//	get index-th lower obstacle value
	value_type lowerObsVal = (*m_spVecOfLowObsValues)[index];

	UG_ASSERT(GetSize(tmpSol) == GetSize(lowerObsVal), "size of tmpSol and size "
			"of lowerObsVal need to be the same");

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
			m_vActiveIndicesLow.push_back(MultiIndex<2>(index, j) );
		}
		else
		{
			//	the 'tmpSol' is valid with respect to the lower constraints
			BlockRef(lastSol[index], j) = BlockRef(tmpSol, j);
		}
	}
}

template <typename TAlgebra>
void
ScalarObstacle<TAlgebra>::
correction_for_upper_obs(vector_type& c, vector_type& lastSol, const size_t index, const value_type& tmpSol)
{
	//	get index-th upper obstacle value
	value_type upperObsVal = (*m_spVecOfUpObsValues)[index];

	UG_ASSERT(GetSize(tmpSol) == GetSize(upperObsVal), "size of tmpSol and size "
			"of upperObsVal need to be the same");

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
			m_vActiveIndicesUp.push_back(MultiIndex<2>(index, j) );
		}
		else
		{
			//	the 'tmpSol' is valid with respect to the upper constraints
			BlockRef(lastSol[index], j) = BlockRef(tmpSol, j);
		}
	}
}

template <typename TAlgebra>
void
ScalarObstacle<TAlgebra>::
correction_for_lower_and_upper_obs(vector_type& c, vector_type& lastSol, const size_t index, const value_type& tmpSol)
{
	//	get index-th lower obstacle value
	value_type upperObsVal = (*m_spVecOfUpObsValues)[index];
	value_type lowerObsVal = (*m_spVecOfLowObsValues)[index];

	UG_ASSERT(GetSize(tmpSol) == GetSize(upperObsVal), "size of tmpSol and size "
			"of upperObsVal need to be the same");
	UG_ASSERT(GetSize(tmpSol) == GetSize(lowerObsVal), "size of tmpSol and size "
			"of lowerObsVal need to be the same");

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
			m_vActiveIndicesUp.push_back(MultiIndex<2>(index, j) );
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
				m_vActiveIndicesLow.push_back(MultiIndex<2>(index, j) );
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

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__SCALAR_OBSTACLE_IMPL__ */
