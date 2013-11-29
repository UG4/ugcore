/*
 * signorini_obstacle.h
 *
 *  Created on: 26.11.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__OBSTACLE_IN_NORMAL_DIR__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__OBSTACLE_IN_NORMAL_DIR__

#include "obstacle_constraint_interface.h"
#include "lib_disc/function_spaces/grid_function.h"

namespace ug{

/// Obstacle Class for Obstacle in normal direction
/**
 *  Scalar obstacle are described by constraints of the form
 *
 * 			u * n <= upObs 		(cf. 'set_upper_obstacle' in 'IObstacleConstraint')
 *
 * 	and
 *
 * 			u * n >= lowObs 	(cf. 'set_lower_obstacle' in 'IObstacleConstraint')
 *
 * 	where u is the solution vector and n the normal vector to the boundary-face/side.
 *
 * 	Here, 'upObs' and 'lowObs' are user-defined functions,
 * 	which need to be of the same size as the function of unknowns u.
 *
 * 	Those obstacle functions can be used in combination with projected preconditioners. They
 * 	should be passed to the preconditioner by 'IProjPreconditioner::set_obstacle_constraint'.
 */
template <typename TDomain, typename TAlgebra>
class ObstacleInNormalDir:
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

	///	Type of grid function
		typedef GridFunction<TDomain, TAlgebra> function_type;

	public:
	/// constructor for an obstacle in normal direction
	///	defined on some subset(s)
		ObstacleInNormalDir(const function_type& u, const char* subsets):
			IObstacleConstraint<TAlgebra>(u, subsets){};

	/// constructor
		ObstacleInNormalDir(): IObstacleConstraint<TAlgebra>(){};

	///	computes the correction for the case that only a lower obstacle is set, i.e. u * n >= g_low
		void correction_for_lower_obs(vector_type& c, vector_type& lastSol, const size_t index);

	///	computes the correction for the case that only an upper obstacle is set, i.e. u * n <= g_up
		void correction_for_upper_obs(vector_type& c, vector_type& lastSol, const size_t index);

	///	computes the correction for the case that a lower and an upper obstacle is set
		void correction_for_lower_and_upper_obs(vector_type& c, vector_type& lastSol, const size_t index);

	///	Destructor
		~ObstacleInNormalDir(){};

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
#include "obstacle_in_normal_dir_impl.h"

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__OBSTACLE_IN_NORMAL_DIR__ */

