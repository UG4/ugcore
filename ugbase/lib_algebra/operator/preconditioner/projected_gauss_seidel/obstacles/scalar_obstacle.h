/*
 * scalar_lower_obstacle.h
 *
 *  Created on: 25.11.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__SCALAR_LOWER_OBSTACLE__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__SCALAR_LOWER_OBSTACLE__

#include "obstacle_constraint_interface.h"
#include "lib_disc/function_spaces/grid_function.h"

using namespace std;

namespace ug{

/// Scalar Lower Obstacles
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
template <typename TDomain, typename TAlgebra>
class ScalarLowerObstacle:
	public IObstacleConstraint<TDomain,TAlgebra>
{
	public:
	///	Base class type
		typedef IObstacleConstraint<TDomain,TAlgebra> base_type;

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
	/// constructor for a scalar obstacle
		ScalarLowerObstacle(const function_type& u):
			IObstacleConstraint<TDomain,TAlgebra>(u){};

	///	default constructor
		ScalarLowerObstacle():
			IObstacleConstraint<TDomain,TAlgebra>(){};

	///	projects the i-th index of the solution onto the admissible set and adjusts the correction
		void adjust_sol_and_cor(value_type& sol_i, value_type& c_i, bool& dofIsAdmissible,
				const DoFIndex& dof);

		void adjust_defect(vector_type& d);

		void restrict_obs_values();

	///	Destructor
		~ScalarLowerObstacle(){};

	private:
	///	store the dofs, which satisfy the constraints with equality
		using base_type::m_vActiveDofs;

	///	map storing the obstacle values for every obstacle dof (key)
		using base_type::m_mObstacleValues;
};

template <typename TDomain, typename TAlgebra>
class ScalarUpperObstacle:
	public IObstacleConstraint<TDomain,TAlgebra>
{
	public:
	///	Base class type
		typedef IObstacleConstraint<TDomain,TAlgebra> base_type;

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
	/// constructor for a scalar obstacle
		ScalarUpperObstacle(const function_type& u):
			IObstacleConstraint<TDomain,TAlgebra>(u){};

	///	default constructor
		ScalarUpperObstacle():
			IObstacleConstraint<TDomain,TAlgebra>(){};

	///	projects the i-th index of the solution onto the admissible set and adjusts the correction
		void adjust_sol_and_cor(value_type& sol_i, value_type& c_i, bool& dofIsAdmissible,
				const DoFIndex& dof);

		void adjust_defect(vector_type& d);

		void restrict_obs_values();

	///	Destructor
		~ScalarUpperObstacle(){};

	private:
	///	store the dofs, which satisfy the constraints with equality
		using base_type::m_vActiveDofs;

	///	map storing the obstacle values for every obstacle dof (key)
		using base_type::m_mObstacleValues;
};

} // end namespace ug

// include implementation
#include "scalar_obstacle_impl.h"

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__SCALAR_LOWER_OBSTACLE__ */
