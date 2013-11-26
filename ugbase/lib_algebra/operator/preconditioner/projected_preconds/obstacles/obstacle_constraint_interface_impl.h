/*
 * obstacle_constraint_interface_impl.h
 *
 *  Created on: 26.11.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__OBSTACLE_CONSTRAINT_INTERFACE_IMPL__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__OBSTACLE_CONSTRAINT_INTERFACE_IMPL__

#include "obstacle_constraint_interface.h"

namespace ug{

template <typename TAlgebra>
void
IObstacleConstraint<TAlgebra>::init(const vector_type& u)
{
// 	init vector of obstacle values and init values with zero
	if (!m_bLowerObs)
		m_spVecOfLowObsValues = u.clone_without_values();
	if (!m_bUpperObs)
		m_spVecOfUpObsValues = u.clone_without_values();

//	check, that lower obstacle is <= upper obstacle (for all indices)
	if (m_bLowerObs && m_bUpperObs)
	{
		if ((*m_spVecOfLowObsValues).size() != (*m_spVecOfUpObsValues).size())
			UG_THROW("In IObstacleConstraint::init(u) :Vector of lower obstacle values [size= "
					<<(*m_spVecOfLowObsValues).size()<<"] and "
					" Vector of upper obstacle values [size= "
					<<(*m_spVecOfUpObsValues).size()<<"] sizes have to match!");

		for(size_t i = 0; i < (*m_spVecOfLowObsValues).size(); i++)
		{
			value_type lowerObsVal = (*m_spVecOfLowObsValues)[i];
			value_type upperObsVal = (*m_spVecOfUpObsValues)[i];
			for(size_t j = 0; j < GetSize(lowerObsVal); j++)
			{
				if (BlockRef(lowerObsVal, j) - BlockRef(upperObsVal, j) > 0.0)
					UG_THROW("In IObstacleConstraint::init(u) " <<i<<"-th index and "<<j<<"-th"
						" component of vector of lower obstacle [value= "<<lowerObsVal<<"] needs "
						"to be lower equal the "<<i<<"-th value of vector of upper obstacle "
						"[value= "<<upperObsVal<<"]!");
			}
		}
	}
}

} // end namespace ug

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__OBSTACLE_CONSTRAINT_INTERFACE_IMPL__ */
