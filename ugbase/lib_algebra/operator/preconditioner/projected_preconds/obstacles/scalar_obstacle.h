/*
 * scalar_obstacle.h
 *
 *  Created on: 25.11.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__SCALAR_OBSTACLE__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__SCALAR_OBSTACLE__

#include "obstacle_constraint_interface.h"

namespace ug{

/// Scalar Obstacles
/**
 *  Scalar obstacle are described by constraints of the form
 *
 * 			u <= upObs 		(cf. 'set_upper_obstacle' in 'IObstacleConstraint')
 *
 * 	and
 *
 * 			u >= lowObs 	(cf. 'set_lower_obstacle' in 'IObstacleConstraint')
 *
 * 	where u is the solution vector. Here, 'upObs' and 'lowObs' are user-defined functions,
 * 	which need to be of the same size as the function of unknowns u.
 *
 * 	Those obstacle functions can be used in combination with projected preconditioners. They
 * 	should be passed to the preconditioner by 'IProjPreconditioner::set_obstacle_constraint'.
 */
template <typename TAlgebra>
class ScalarObstacle:
	public IObstacleConstraint<TAlgebra>
{
	public:
	///	Base class type
		typedef IObstacleConstraint<TAlgebra> base_type;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Matrix type
		typedef typename algebra_type::matrix_type matrix_type;

	///	Vector type
		typedef typename algebra_type::vector_type vector_type;

	///	Value type
		typedef typename vector_type::value_type value_type;

	public:
	/// constructor
		ScalarObstacle(): IObstacleConstraint<TAlgebra>(){};

	///	computes the correction for the case that only a lower obstacle is set, i.e. u >= g_low
		void correction_for_lower_obs(vector_type& c, vector_type& lastSol, const size_t index, const value_type& tmpSol);

	///	computes the correction for the case that only an upper obstacle is set, i.e. u <= g_up
		void correction_for_upper_obs(vector_type& c, vector_type& lastSol, const size_t index, const value_type& tmpSol);

	///	computes the correction for the case that a lower and an upper obstacle is set
		void correction_for_lower_and_upper_obs(vector_type& c, vector_type& lastSol, const size_t index, const value_type& tmpSol);

	///	Destructor
		~ScalarObstacle(){};

	private:
	///	pointer to vector of active Indices (for lower and upper constraint)
		using base_type::m_spLowerActiveInd;
		using base_type::m_spUpperActiveInd;

	///	vector of obstacle values (for lower and upper constraint)
		using base_type::m_spVecOfLowObsValues;
		using base_type::m_spVecOfUpObsValues;
};

} // end namespace ug

// include implementation
#include "scalar_obstacle_impl.h"

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__SCALAR_OBSTACLE__ */
