/*
 * Copyright (c) 2018:  G-CSC, Goethe University Frankfurt
 * Author: Martin Stepniewski
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

#ifndef __H__UG__LIB_ALGEBRA__POWER_METHOD_H__
#define __H__UG__LIB_ALGEBRA__POWER_METHOD_H__

#include <vector>
#include <string>
#include "additional_math.h"

#include "lib_algebra/operator/interface/linear_operator.h"
#include "lib_algebra/operator/interface/linear_operator_inverse.h"
#include "lib_algebra/operator/interface/preconditioner.h"
#include "lib_algebra/operator/interface/matrix_operator.h"

#include "lib_algebra/algebra_common/vector_util.h"
#include  "common/profiler/profiler.h"


namespace ug{


/**
 * Power Method Eigensolver
 *
 * The Power Method solves generalized eigenvalue problems of the form A v = lambda B v
 * calculating the largest/smallest (in terms of abs(lambda)) eigenvalues of the problem
 * for sparse matrices A and B which e.g. emerge from FE or FV discretizations.
 */
template<typename TAlgebra>
class PowerMethod
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;
		typedef typename TAlgebra::matrix_type matrix_type;
		typedef typename vector_type::value_type vector_value_type;

	private:
	//	Stiffness matrix
		SmartPtr<ILinearOperator<vector_type> > m_spLinOpA;

	//	Mass matrix
		typedef typename IPreconditioner<TAlgebra>::matrix_operator_type matrix_operator_type;
		SmartPtr<ILinearOperator<vector_type> > m_spLinOpB;
		SmartPtr<matrix_operator_type> m_spMatOpA;
		SmartPtr<matrix_operator_type> m_spMatOpB;
		//matrix_operator_type *m_pB;

	//	Eigenvector and -value
		SmartPtr<vector_type> m_spEigenvector;
		SmartPtr<vector_type> m_spEigenvectorOld;
		double m_dMaxEigenvalue;
		double m_dMinEigenvalue;

	//	Solver
		SmartPtr<ILinearOperatorInverse<vector_type> > m_spSolver;

	//	Precision
		size_t m_maxIterations;
		double m_dPrecision;

	//	Residual
		SmartPtr<vector_type> m_spResidual;

	//	Dirichlet node indexing
		std::vector<bool> m_vbDirichlet;
		int m_numDirichletRows;

	public:
	//	Number of iterations
		size_t m_iteration;

		PowerMethod()
		{
			UG_LOG("Initializing PowerMethod." << std::endl);
			m_spLinOpA = SPNULL;
			m_spLinOpB = SPNULL;
			m_spMatOpA = SPNULL;
			m_spMatOpB = SPNULL;
			m_maxIterations = 1000;
			m_dPrecision = 1e-8;
			m_iteration = 0;
			m_dMaxEigenvalue = 0.0;
			m_dMinEigenvalue = 0.0;
			m_numDirichletRows = 0;
		}

		void set_solver(SmartPtr<ILinearOperatorInverse<vector_type> > solver)
		{
			m_spSolver = solver;
		}

		void set_linear_operator_A(SmartPtr<ILinearOperator<vector_type> > loA)
		{
			m_spLinOpA = loA;

		// 	get dirichlet nodes
			m_spMatOpA = m_spLinOpA.template cast_dynamic<MatrixOperator<matrix_type, vector_type> >();
			matrix_type& A = m_spMatOpA->get_matrix();
			m_vbDirichlet.resize(A.num_rows());

			for(size_t i=0; i<A.num_rows(); i++)
			{
				m_vbDirichlet[i] = A.is_isolated(i);

				if(A.is_isolated(i) == 0)
					m_numDirichletRows++;
			}
		}

		void set_linear_operator_B(SmartPtr<ILinearOperator<vector_type> > loB)
		{
			m_spLinOpB = loB;
			m_spMatOpB = m_spLinOpB.template cast_dynamic<MatrixOperator<matrix_type, vector_type> >();
		}

		void set_start_vector(SmartPtr<vector_type> vec)
		{
			m_spEigenvector = vec;
		}

		void set_max_iterations(size_t maxIterations)
		{
			m_maxIterations = maxIterations;
		}

		void set_precision(double precision)
		{
			m_dPrecision = precision;
		}

		int calculate_max_eigenvalue()
		{
			PROFILE_FUNC_GROUP("PowerMethod");

			UG_COND_THROW(m_spSolver == SPNULL && m_spLinOpB != SPNULL, "PowerMethod::calculate_max_eigenvalue(): Solver not set, please specify.");
			if(m_spLinOpB != SPNULL)
				m_spSolver->init(m_spLinOpB);

			m_spResidual = create_approximation_vector();

			for(m_iteration = 0; m_iteration < m_maxIterations; ++m_iteration)
			{
				m_spEigenvectorOld = m_spEigenvector->clone();

			//	residual = A v
				UG_COND_THROW(m_spLinOpA == SPNULL, "PowerMethod::calculate_max_eigenvalue(): Linear operator A not set, please specify.");
				m_spLinOpA->apply(*m_spResidual, *m_spEigenvector);

			//	v = B^-1 A v
				if(m_spLinOpB != SPNULL)
					m_spSolver->apply(*m_spEigenvector, *m_spResidual);
			//	in case B is not set: v = A v
				else
					m_spEigenvector = m_spResidual->clone();

			//	v = v / ||v||_B
				normalize_approximations();

			//	in case B is not set, restore storage type of m_spEigenvector to PST_CONSISTENT
			//	changed by norm calculation in normalize_approximations()
				if(m_spLinOpB == SPNULL)
				{
					#ifdef UG_PARALLEL
							m_spEigenvector->change_storage_type(PST_CONSISTENT);
					#endif
				}

			//	residual = v - v_old
				VecScaleAdd(*m_spResidual, 1.0, *m_spEigenvectorOld, -1.0, *m_spEigenvector);

				if(m_spResidual->norm() <= m_dPrecision)
				{
					UG_LOG("PowerMethod::calculate_max_eigenvalue() converged after " << m_iteration << " iterations." << std::endl);
					break;
				}

				if(m_iteration == m_maxIterations-1)
					UG_LOG("PowerMethod::calculate_max_eigenvalue() reached precision of " << m_spResidual->norm() << " after " << m_maxIterations << " iterations." << std::endl);
			}

		//	lambda = <v,Av>
			m_spLinOpA->apply(*m_spResidual, *m_spEigenvector);
			m_dMaxEigenvalue = m_spEigenvector->dotprod(*m_spResidual);

			return true;
		}

		int calculate_min_eigenvalue()
		{
			PROFILE_FUNC_GROUP("PowerMethod");

			UG_COND_THROW(m_spLinOpA == SPNULL, "PowerMethod::calculate_min_eigenvalue(): Linear operator A not set, please specify.");
			UG_COND_THROW(m_spSolver == SPNULL, "PowerMethod::calculate_min_eigenvalue(): Solver not set, please specify.");

			m_spResidual = create_approximation_vector();

			m_spSolver->init(m_spLinOpA);

			for(m_iteration = 0; m_iteration < m_maxIterations; ++m_iteration)
			{
				m_spEigenvectorOld = m_spEigenvector->clone();

			//	residual = B v
				if(m_spLinOpB != SPNULL)
					m_spLinOpB->apply(*m_spResidual, *m_spEigenvector);
			//	in case B is not set
				else
				{
					m_spResidual = m_spEigenvector->clone();

					#ifdef UG_PARALLEL
						m_spResidual->change_storage_type(PST_ADDITIVE);
					#endif
				}

			//	v = A^-1 B v or v = A^-1 v
				m_spSolver->apply(*m_spEigenvector, *m_spResidual);

			//	v = v / ||v||_B
				normalize_approximations();

			//	in case B is not set, restore storage type of m_spEigenvector to PST_CONSISTENT
			//	changed by norm calculation in normalize_approximations()
				if(m_spLinOpB == SPNULL)
				{
					#ifdef UG_PARALLEL
							m_spEigenvector->change_storage_type(PST_CONSISTENT);
					#endif
				}

			//	residual = v - v_old
				VecScaleAdd(*m_spResidual, 1.0, *m_spEigenvectorOld, -1.0, *m_spEigenvector);

				if(m_spResidual->norm() <= m_dPrecision)
				{
					UG_LOG("PowerMethod::calculate_min_eigenvalue() converged after " << m_iteration << " iterations." << std::endl);
					break;
				}

				if(m_iteration == m_maxIterations-1)
					UG_LOG("PowerMethod::calculate_min_eigenvalue() reached precision of " << m_spResidual->norm() << " after " << m_maxIterations << " iterations." << std::endl);
			}

		//	lambda = <v,Av>
			m_spLinOpA->apply(*m_spResidual, *m_spEigenvector);
			m_dMinEigenvalue = m_spEigenvector->dotprod(*m_spResidual);

			return true;
		}

		double get_max_eigenvalue()
		{
			return m_dMaxEigenvalue;
		}

		double get_min_eigenvalue()
		{
			return m_dMinEigenvalue;
		}

		size_t get_iterations()
		{
			return m_iteration;
		}

		void print_matrix_A()
		{
			PrintMatrix(m_spMatOpA->get_matrix(), "A");
		}

		void print_matrix_B()
		{
			PrintMatrix(m_spMatOpB->get_matrix(), "B");
		}

		void print_eigenvector()
		{
			m_spEigenvector->print();
		}

	private:

		SmartPtr<vector_type> create_approximation_vector()
		{
			SmartPtr<vector_type> t = m_spEigenvector->clone_without_values();
			t->set(0.0);
			return t;
		}

		double B_norm(vector_type &x)
		{
			if(m_spMatOpB != SPNULL)
				return EnergyNorm(x, *m_spMatOpB);
			else
				return x.norm();
		}

		void normalize_approximations()
		{
			*m_spEigenvector *= 1/ B_norm(*m_spEigenvector);
		}
};

} // namespace ug


#endif // __H__UG__LIB_ALGEBRA__POWER_METHOD_H__
