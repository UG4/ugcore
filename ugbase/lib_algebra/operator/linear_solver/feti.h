/*
 * feti.h
 *
 *  Created on: 11.11.2010
 *      Author: iheppner, avogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__FETI__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__FETI__

#include <iostream>
#include <sstream>
#include "lib_algebra/operator/operator_inverse_interface.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

template <typename TAlgebra>
class FETISolver : public IMatrixOperatorInverse<	typename TAlgebra::vector_type,
													typename TAlgebra::vector_type,
													typename TAlgebra::matrix_type>
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
	///	constructor
		FETISolver() :
			m_A(NULL), m_pSequentialSolver(NULL), m_pConvCheck(NULL)
		{}

	///	name of solver
		virtual const char* name() const {return "FETI Solver";}

	///	sets a convergence check
		void set_convergence_check(IConvergenceCheck& convCheck)
		{
			m_pConvCheck = &convCheck;
			m_pConvCheck->set_offset(3);
		}

	/// returns the convergence check
		IConvergenceCheck* get_convergence_check() {return m_pConvCheck;}

	///	sets a sequential solver
		void set_solver(IMatrixOperatorInverse<vector_type, vector_type, matrix_type>& seqSolver)
		{
			m_pSequentialSolver = &seqSolver;
		}

	///	initializes the solver for operator A
		virtual bool init(IMatrixOperator<vector_type, vector_type, matrix_type>& A)
		{
		//	save current operator
			m_A = &A;

		//	save matrix to invert
			m_pMatrix = m_A->get_matrix();

		//	init sequential solver
			if(m_pSequentialSolver != NULL)
				if(!m_pSequentialSolver->init(m_A))
				{
					UG_LOG("ERROR in 'FETISolver::prepare': Cannot init "
							"Sequential Linear Solver for Operator A.\n");return false;
				}

		//	we're done
			return true;
		}

	///	solves the system and returns the last defect of iteration in rhs
		virtual bool apply_return_defect(vector_type& x, vector_type& b)
		{
		//	check that matrix has been set
			if(m_A == NULL)
			{
				UG_LOG("ERROR: In 'FETISolver::apply': Matrix A not set.\n");
				return false;
			}

		//	check that sequential solver has been set
			if(m_pSequentialSolver != NULL)
			{
				UG_LOG("ERROR: In 'FETISolver::apply': No Sequential Solver set.\n");
				return false;
			}

		//	check that convergence check is set
			if(m_pConvCheck == NULL)
			{
				UG_LOG("ERROR: In 'FETISolver::apply':"
						" Convergence check not set.\n");
				return false;
			}

		//	Check parallel storage type of vectors
			#ifdef UG_PARALLEL
			if(!b.has_storage_type(PST_ADDITIVE) || !x.has_storage_type(PST_CONSISTENT))
			{
				UG_LOG("ERROR: In 'FETISolver::apply': "
						"Inadequate storage format of Vectors.\n");
				return false;
			}
			#endif

		// 	Rename b as d (for convenience)
			vector_type& d = b;

		// 	Build defect:  d := b - A*x
			if(!m_A->apply_sub(d, x))
			{
				UG_LOG("ERROR in 'LinearOperatorInverse::apply': "
						"Unable to build defect. Aborting.\n");
				return false;
			}

		// 	Create correction
		// 	todo: 	it would be sufficient to only copy the pattern (and parallel constructor)
		//			without initializing the values
			vector_type c; c.create(x.size()); c = x;

		//  Prepare convergence check
			prepare_conv_check();

		//	Compute norm of first defect (in parallel)
			m_pConvCheck->start(d);

		// 	Iteration loop
			while(!m_pConvCheck->iteration_ended())
			{

				// TODO: Implement Solver

			// 	Add (parallel) correction to solution
				x += c;

			// 	compute new defect (in parallel)
				m_pConvCheck->update(d);
			}

		//	Post Output
			if(!m_pConvCheck->post())
			{
				UG_LOG("ERROR in 'LinearSolver::apply': "
						"post-convergence-check signaled failure. Aborting.\n");
				return false;
			}

		//	we're done
			return true;
		}

	///	solves the system
		virtual bool apply(vector_type& x, const vector_type& b)
		{
		//	copy defect
			vector_type d; d.resize(b.size());
			d = b;

		//	solve on copy of defect
			return apply_return_defect(x, d);
		}

		// destructor
		virtual ~LinearSolver() {};

	protected:
	//	Prepare the convergence check
		void prepare_conv_check()
		{
			m_pConvCheck->set_name(name());
			m_pConvCheck->set_symbol('%');
			m_pConvCheck->set_name(name());

			if(m_pSequentialSolver != NULL)
			{
				stringstream ss; ss <<  " (Seq. Solver: " << m_pSequentialSolver->name() << ")";
				m_pConvCheck->set_info(ss.str());
			}
			else
			{
				m_pConvCheck->set_info(" (No Seq. Solver) ");
			}
		}

	protected:
	// 	Operator that is inverted by this Inverse Operator
		ILinearOperator<vector_type,vector_type>* m_A;

	// 	Parallel Matrix to invert
		matrix_type* m_pMatrix;

	// 	Linear Solver to invert the process local problems
		IMatrixOperatorInverse<vector_type,vector_type, matrix_type>* m_pSequentialSolver;

	// 	Convergence Check
		IConvergenceCheck* m_pConvCheck;
};

} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__FETI__ */
