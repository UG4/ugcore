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

	///	Domain type
		typedef TDomain domain_type;

	///	World dimension
		static const int dim = domain_type::dim;

	///	Type of position coordinates (e.g. position_type)
		typedef typename domain_type::position_type position_type;

	public:
	/// constructor for an obstacle in normal direction
	///	defined on some subset(s)
		ObstacleInNormalDir(const function_type& u):
			IObstacleConstraint<TDomain,TAlgebra>(u){
			m_spDD = u.dof_distribution();
			m_spDomain = u.domain();

			UG_LOG("In ObstacleInNormalDir::constructor u hat "<<u.size()<<"Eintraege \n");
			UG_LOG("\n");
		};

	/// constructor
		ObstacleInNormalDir(): IObstacleConstraint<TDomain,TAlgebra>(){};

	///	preprocess is useful to attach the normals for every obstacle DoF
		void preprocess();

	///	projects the i-th index of the solution onto the admissible set and adjusts the correction
		void adjust_sol_and_cor(value_type& sol_i, value_type& c_i, bool& dofIsAdmissible,
				const DoFIndex& dof);

		void adjust_defect(vector_type& d);

		void restrict_obs_values();

	///	Destructor
		~ObstacleInNormalDir(){};

	private:
		void transform_eulerian_coord_sys(MathVector<dim> transformedONB[],
				const MathVector<dim>& firstTransformedBaseVec);

		template <typename TElem, typename TIterator>
		void adjust_sol_and_cor_elem(TIterator iterBegin, TIterator iterEnd, value_type& sol_i,
				value_type& c_i, bool& dofIsAdmissible, const DoFIndex& dof);

	private:
	///	store the dofs, which satisfy the constraints with equality
		using base_type::m_vActiveDofs;

	///	map storing the obstacle values for every obstacle dof (key)
		using base_type::m_mObstacleValues;

	///	stores the subset-indices of the obstacle subsets
		using base_type::m_vObsSubsets;

	///	pointer to the DofDistribution on the whole domain
		ConstSmartPtr<DoFDistribution> m_spDD;

	///	pointer to the domain
		ConstSmartPtr<TDomain> m_spDomain;

	///	struct to store data for a specific obstacle DoF
		struct obsDoFData{
			number obsVal;
			MathVector<dim> transformedONB[dim];
		};
};

} // end namespace ug

// include implementation
#include "obstacle_in_normal_dir_impl.h"

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__OBSTACLE_IN_NORMAL_DIR__ */

