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
		using base_type = IObstacleConstraint<TDomain,TAlgebra>;

	///	Algebra type
		using algebra_type = TAlgebra;

	///	Matrix type
		using matrix_type = typename algebra_type::matrix_type;

	///	Vector type
		using vector_type = typename algebra_type::vector_type;

	///	Value type
		using value_type = typename vector_type::value_type;

	///	Type of grid function
		using function_type = GridFunction<TDomain, TAlgebra>;

	public:
	/// constructor for a scalar obstacle
		ScalarLowerObstacle(const function_type& u):
			IObstacleConstraint<TDomain,TAlgebra>(u){};

	///	default constructor
		ScalarLowerObstacle():
			IObstacleConstraint<TDomain,TAlgebra>(){};

	///	projects the i-th index of the solution onto the admissible set and adjusts the correction
		void adjust_sol_and_cor(value_type& sol_i, value_type& c_i, bool& dofIsActive,
				const DoFIndex& dof);

		void adjust_defect_to_constraint(vector_type& d);

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
		using base_type = IObstacleConstraint<TDomain,TAlgebra>;

	///	Algebra type
		using algebra_type = TAlgebra;

	///	Matrix type
		using matrix_type = typename algebra_type::matrix_type;

	///	Vector type
		using vector_type = typename algebra_type::vector_type;

	///	Value type
		using value_type = typename vector_type::value_type;

	///	Type of grid function
		using function_type = GridFunction<TDomain, TAlgebra>;

	public:
	/// constructor for a scalar obstacle
		ScalarUpperObstacle(const function_type& u):
			IObstacleConstraint<TDomain,TAlgebra>(u){};

	///	default constructor
		ScalarUpperObstacle():
			IObstacleConstraint<TDomain,TAlgebra>(){};

	///	projects the i-th index of the solution onto the admissible set and adjusts the correction
		void adjust_sol_and_cor(value_type& sol_i, value_type& c_i, bool& dofIsActive,
				const DoFIndex& dof);

		void adjust_defect_to_constraint(vector_type& d);

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
