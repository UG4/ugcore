/*
 * feti.h
 *
 *  Created on: 11.11.2010
 *      Author: iheppner, avogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__FETI__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__FETI__

namespace ug{

#ifdef UG_PARALLEL

#include <iostream>
#include <sstream>
#include <string>
#include "lib_algebra/operator/operator_inverse_interface.h"
#include "lib_algebra/parallelization/parallelization.h"
#include "lib_algebra/operator/debug_writer.h"

/// operator implementation of the local Schur complement
/**
 * This operator is the application of the local Schur complement. The underlying
 * matrix must have at least two layouts. The first layout will be used as
 * inner parts, all other layouts are used to identify the boundary and the
 * schur complement is build w.r.t. to these variables.
 */
template <typename TAlgebra>
class LocalSchurComplement
	: public ILinearOperator<	typename TAlgebra::vector_type,
	  	  	  	  	  	  	  	typename TAlgebra::vector_type>
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
		LocalSchurComplement()
		{}

	///	name of solver
		virtual const char* name() const {return "Local Schur Complement Solver";}

	///	sets a sequential Dirichlet solver
	/**
	 * This method sets the Dirichlet Solver that is used to invert the
	 * inner matrix \f$A_{II}\F$
	 */
		void set_dirichlet_solver(ILinearOperatorInverse<vector_type, vector_type>& dirichletSolver)
		{
		//	remember the Dirichlet Solver
			m_pDirichletSolver = &dirichletSolver;
		}

	///	set debug output
		void set_debug(IDebugWriter<algebra_type>* debugWriter)
		{
			m_pDebugWriter = debugWriter;
		}

	///	set original matrix from which the local Schur complement is constructed
	/**
	 * Using this method, the original matrix A is set. Given the matrix in the
	 * form \f$ A = \sum\limits_{p=1}^{N} A^{(p)}\f$ in additive form, with
	 * \f{align*}
	 * A^{(p)}
	 * \begin{pmatrix}
	 * A_{II}^{(p)} & A_{I \Delta}^{(p)} & A_{I \Pi}^{(p)} \\
	 * A_{\Delta I}^{(p)} & A_{\Delta \Delta}^{(p)} & A_{\Delta \Pi}^{(p)} \\
	 * A_{\Pi I}^{(p)} & A_{\Pi \Delta}^{(p)} & A_{\Pi \Pi}^{(p)}
	 * \end{pmatrix}
	 * \f}
	 *
	 * the local schur complement is the processwise application of the
	 * operator
	 *
	 * \f{align*}
	 * S_{\Delta \Delta} = A_{\Delta \Delta} -
	 * \begin{pmatrix} A_{I \Delta}^T & A_{\Pi \Delta}^T \end{pmatrix}
	 * \begin{pmatrix} A_{II} & A_{I \Pi} \\ A_{\Pi I} & A_{\Pi \Pi} \end{pmatrix}^{-1}
	 * \begin{pmatrix} A_{I \Delta} \\ A_{\Pi \Delta} \end{pmatrix}	 *
	 * \f}
	 */
		void set_matrix(IMatrixOperator<vector_type, vector_type, matrix_type>& A)
		{
		//	save current operator
			m_pOperator = &A;
		}

	/// implementation of the operator for the solution dependent initialization.
		bool init(const vector_type& u) {return init();}

	///	initializes the solver for operator A
	/**
	 * This method must be called, before the apply() method can be invoked.
	 * It has to be called each time, when the matrix has been replaced. A deep
	 * copy of the matrix is then constructed and in this copy the rows belonging
	 * to the \f$\Delta\f$ and \f$\Pi\f$ unknowns are set to identity rows. This
	 * matrix is used in the solution of the local dirichlet problem.
	 */
		virtual bool init()
		{
		//	check that operator has been set
			if(m_pOperator == NULL)
			{
				UG_LOG("ERROR in 'LocalSchurComplement::init': No Operator A"
						" set.\n");
				return false;
			}

		//	save matrix to invert
			m_pMatrix = &m_pOperator->get_matrix();

		//	check that matrix has enough decomposition level
			if(m_pMatrix->num_layouts() < 2)
			{
				UG_LOG("ERROR in 'LocalSchurComplement::init': The Operator must"
						" have at least two layouts, but the current Operator has"
						" only " << m_pMatrix->num_layouts() << "\n");
				return false;
			}

		//	get matrix from dirichlet operator
			m_pDirichletMatrix = &m_DirichletOperator.get_matrix();

		//	Copy Matrix for Dirichlet Problem
			*m_pDirichletMatrix = *m_pMatrix;

		//	Set Dirichlet values on Delta and Pi
			for(size_t i = 1; i < m_pDirichletMatrix->num_layouts(); ++i)
			{
				MatSetDirichletOnLayout(m_pDirichletMatrix, m_pDirichletMatrix->get_slave_layout(i));
				MatSetDirichletOnLayout(m_pDirichletMatrix, m_pDirichletMatrix->get_master_layout(i));
			}

		//	init sequential solver for Dirichlet problem
			if(m_pDirichletSolver != NULL)
				if(!m_pDirichletSolver->init(m_DirichletOperator))
				{
					UG_LOG("ERROR in 'LocalSchurComplement::init': Cannot init "
							"Sequential Dirichlet Solver for Operator A.\n");return false;
				}

		//	Debug output of matrices
			if(m_pDebugWriter != NULL)
			{
				m_pDebugWriter->write_matrix(m_DirichletOperator.get_matrix(),
				                             "FetiDirichletMatrix");
				m_pDebugWriter->write_matrix(m_pOperator->get_matrix(),
				                             "FetiNeumannMatrix");
			}

		//	we're done
			return true;
		}

	///	solves the system and returns the last defect of iteration in rhs
		virtual bool apply(vector_type& f, const vector_type& u)
		{
		//	check that matrix has been set
			if(m_pOperator == NULL)
			{
				UG_LOG("ERROR: In 'LocalSchurComplement::apply': "
								"Matrix A not set.\n");
				return false;
			}

		//	check dirichlet solver
			if(m_pDirichletSolver == NULL)
			{
				UG_LOG("ERROR: In 'LocalSchurComplement::apply':"
								" No sequential Dirichlet Solver set.\n");
				return false;
			}

		//	Check parallel storage type of vectors
			if(!f.has_storage_type(PST_ADDITIVE) || !u.has_storage_type(PST_CONSISTENT))
			{
				UG_LOG("ERROR: In 'LocalSchurComplement::apply': "
								"Inadequate storage format of Vectors.\n");
				return false;
			}

		//	Help vector
		//	\todo: it would be sufficient to copy only the layouts without copying the values
			vector_type uTmp; uTmp.create(u.size()); uTmp = u;

		//	1. Set values to zero on I and Pi
			// a) Reset all values
			uTmp.set(0.0);

			// b) Copy values on \Delta
			VecScaleAddOnLayout(&uTmp, &u, 1.0, u.get_slave_layout(2));
			VecScaleAddOnLayout(&uTmp, &u, 1.0, u.get_master_layout(2));

		//	2. Compute rhs f_{I} = A_{I \Delta} u_{\Delta}
			if(!m_DirichletOperator.apply(f, uTmp))
			{
				UG_LOG_ALL_PROCS("ERROR in 'LocalSchurComplement::apply': "
								 "Could not compute Rhs for Dirichlet problem on "
								 "Proc " << pcl::GetProcRank() << ".\n");
				return false;
			}
			// todo: set values to zero on \Delta

		//	3. Invert on inner unknowns u_{I} = A_{II}^{-1} f_{I}
			// a) use the inner-FETI-block layouts
			uTmp.use_layout(0);
			f.use_layout(0);
			m_pDirichletMatrix->use_layout(0);

			// b) invoke Dirichlet solver
			if(!m_pDirichletSolver->apply_return_defect(uTmp, f))
			{
				UG_LOG_ALL_PROCS("ERROR in 'LocalSchurComplement::apply': "
								 "Could not solve Dirichlet problem on Proc "
									<< pcl::GetProcRank() << ".\n");
				return false;
			}

		//	4. Compute result vector
			// a) Scale u_{I} by -1
			uTmp *= -1.0;

			// b) Add u_{\Delta} on \Delta
			VecScaleAddOnLayout(&uTmp, &u, 1.0, u.get_slave_layout(2));
			VecScaleAddOnLayout(&uTmp, &u, 1.0, u.get_master_layout(2));

			// c) Multiply with full matrix
			if(!m_pOperator->apply(f, uTmp))
			{
				UG_LOG_ALL_PROCS("ERROR in 'LocalSchurComplement::apply': "
								 "Could not apply full matrix on "
								 "Proc " << pcl::GetProcRank() << ".\n");
				return false;
			}

		//	5. Reset all values for I, \Pi
			//\todo: implement

		//	we're done
			return true;
		}

	///	solves the system
		virtual bool apply_sub(vector_type& f, const vector_type& u)
		{
		//	create new rhs
			vector_type d; d.resize(f.size());

		//	solve
			if(!apply(d, u)) return false;

		//	subtract from vector
			f -= d;

		//	we're done
			return true;
		}

		// destructor
		virtual ~LocalSchurComplement() {};

	protected:
		bool write_debug(const vector_type& vec, const char* filename)
		{
		//	if no debug writer set, we're done
			if(m_pDebugWriter == NULL) return true;

		//	write
			return m_pDebugWriter->write_vector(vec, filename);
		}

	protected:
	// 	Operator that is inverted by this Inverse Operator
		IMatrixOperator<vector_type,vector_type,matrix_type>* m_pOperator;

	// 	Parallel Matrix
		matrix_type* m_pMatrix;

	//	Copy of matrix
		PureMatrixOperator<vector_type, vector_type, matrix_type> m_DirichletOperator;

	// 	Parallel Dirichlet Matrix
		matrix_type* m_pDirichletMatrix;

	// 	Linear Solver to invert the local Dirichlet problems
		ILinearOperatorInverse<vector_type,vector_type>* m_pDirichletSolver;

	//	Debug Writer
		IDebugWriter<algebra_type>* m_pDebugWriter;
};

#endif /* UG_PARALLEL */

} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__FETI__ */
