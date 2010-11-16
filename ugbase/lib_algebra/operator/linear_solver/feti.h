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
#include "lib_algebra/operator/operator_inverse_interface.h"
#include "lib_algebra/parallelization/parallelization.h"


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
			m_theta(1.0), m_A(NULL), m_pNeumannSolver(NULL),
			m_pDirichletSolver(NULL), m_pConvCheck(NULL)
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

	///	sets a sequential neumann solver
		void set_neumann_solver(ILinearOperatorInverse<vector_type, vector_type>& neumannSolver)
		{
			m_pNeumannSolver = &neumannSolver;
		}

	///	sets a sequential dirichlet solver
		void set_dirichlet_solver(ILinearOperatorInverse<vector_type, vector_type>& dirichletSolver)
		{
			m_pDirichletSolver = &dirichletSolver;
		}

	///	sets damping factor
		void set_theta(number theta)
		{
			m_theta = theta;
		}

	///	initializes the solver for operator A
		virtual bool init(IMatrixOperator<vector_type, vector_type, matrix_type>& A)
		{
		//	save current operator
			m_A = &A;

		//	save matrix to invert
			m_pMatrix = &m_A->get_matrix();

		//	Copy Matrix for Dirichlet Problem
			matrix_type& dirMat = m_DirichletOperator.get_matrix();
			dirMat = *m_pMatrix;

		//	Set dirichlet values
			set_dirichlet_rows_on_gamma(dirMat);

		//	init sequential solver
			if(m_pNeumannSolver != NULL)
				if(!m_pNeumannSolver->init(*m_A))
				{
					UG_LOG("ERROR in 'FETISolver::prepare': Cannot init "
							"Sequential Neumann Solver for Operator A.\n");return false;
				}

		//	init sequential solver
			if(m_pDirichletSolver != NULL)
				if(!m_pDirichletSolver->init(m_DirichletOperator))
				{
					UG_LOG("ERROR in 'FETISolver::prepare': Cannot init "
							"Sequential Dirichlet Solver for Operator A.\n");return false;
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
			if(m_pNeumannSolver == NULL)
			{
				UG_LOG("ERROR: In 'FETISolver::apply': No Sequential Neumann Solver set.\n");
				return false;
			}
			if(m_pDirichletSolver == NULL)
			{
				UG_LOG("ERROR: In 'FETISolver::apply': No Sequential Dirichlet Solver set.\n");
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
			if(!b.has_storage_type(PST_ADDITIVE) || !x.has_storage_type(PST_CONSISTENT))
			{
				UG_LOG("ERROR: In 'FETISolver::apply': "
						"Inadequate storage format of Vectors.\n");
				return false;
			}

		//	prepare layouts and communicators used for inner and intra Feti-Partition communication
			prepare_layouts_and_communicators(x);

		// 	todo: 	it would be sufficient to only copy the pattern (and parallel constructor)
		//			without initializing the values

		//	Normal Flux at Gamma - Boundary (i.e. inner boundary of feti method)
			vector_type lambda; lambda.create(x.size()); lambda = x;

		//	Flux Correction at Gamma - Boundary
			vector_type eta; eta.create(x.size()); eta = x;

		//	Modified right - hand side for neumann problem
			vector_type ModRhs; ModRhs.create(b.size()); ModRhs = b;
			vector_type ModRhsCopy; ModRhsCopy.create(b.size()); ModRhsCopy = b;

		//  Prepare convergence check
			prepare_conv_check();

		//	flag iff first iterate
			bool first = true;

		//	set lambda to zero
			lambda.set(0.0);

		// 	Iteration loop
			while(true)
			{
			//	Copy right-hand side
				ModRhs = b;

			//	Compute new rhs for neumann problem using lambda on Gamma
				add_flux_to_rhs(ModRhs, lambda);

			//	set local communication
				set_layouts_for_local_feti_partition_communication(x);
				set_layouts_for_local_feti_partition_communication(ModRhs);

			//	Solve neumann problem on feti-partition
				if(!m_pNeumannSolver->apply_return_defect(x, ModRhs))
				{
					UG_LOG("ERROR in 'FETI-Solver::apply': "
							"Could not solve neumann problem on proc " << pcl::GetProcRank() << ".\n");
					return false;
				}
				// todo: Check that all processes solved the problem

			//	set global communication
				set_layouts_for_inter_feti_partition_communication(x);
				set_layouts_for_inter_feti_partition_communication(ModRhs);

			//	Set dirichlet values for Rhs, zero else
				copy_dirichlet_values_and_zero(ModRhs, x);

			//	Compute difference of solution on Gamma
				compute_difference_on_gamma(ModRhs);

				ModRhsCopy = ModRhs;
			//	Compute norm of difference on Gamma
				if(first)
				{
				//	Compute first defect
					m_pConvCheck->start(ModRhsCopy);
					first = false;
				}
				else
				{
				// 	compute new defect (in parallel)
					m_pConvCheck->update(ModRhsCopy);
				}

			//	check if iteration ended
				if(m_pConvCheck->iteration_ended())
					break;

			//	set local communication
				set_layouts_for_local_feti_partition_communication(x);
				set_layouts_for_local_feti_partition_communication(ModRhs);

			//	Solve dirichlet problem on feti-partition
				if(!m_pDirichletSolver->apply_return_defect(x, ModRhs))
				{
					UG_LOG("ERROR in 'FETI-Solver::apply': "
							"Could not solve neumann problem on proc " << pcl::GetProcRank() << ".\n");
					return false;
				}
				// todo: Check that all processes solved the problem

			//	set global communication
				set_layouts_for_inter_feti_partition_communication(x);
				set_layouts_for_inter_feti_partition_communication(ModRhs);

			//	Compute update for lambda: eta = A*x
				m_pMatrix->apply(eta, x);

			//	sum up over all processes
				eta.set_storage_type(PST_ADDITIVE);
				eta.change_storage_type(PST_CONSISTENT);

			// 	Update lambda
				VecScaleAdd(lambda, 1.0, lambda, -m_theta, eta);
			}

		//	Post Output
			if(!m_pConvCheck->post())
			{
				UG_LOG("ERROR in 'FETI-Solver::apply': "
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
		virtual ~FETISolver() {};

	protected:
	//	Prepare the convergence check
		void prepare_conv_check()
		{
			m_pConvCheck->set_name(name());
			m_pConvCheck->set_symbol('%');
			m_pConvCheck->set_name(name());

			if(m_pNeumannSolver != NULL && m_pDirichletSolver != NULL)
			{
				stringstream ss; ss <<  " (Seq. Solver: " << m_pNeumannSolver->name() << ","
									<< m_pDirichletSolver->name() << ")";
				m_pConvCheck->set_info(ss.str());
			}
			else
			{
				m_pConvCheck->set_info(" (No Seq. Solver) ");
			}
		}

	protected:
	//	add lambda on the Gamma Boundary
		void add_flux_to_rhs(vector_type& ModRhs, const vector_type& lambda)
		{
			number scale = 1.0;
			if(pcl::GetProcRank() != 1)
				scale = -1.0;

			VecScaleAddOnLayoutWithoutCommunication(&ModRhs, &lambda, scale, m_InterSlaveIndexLayout);
			VecScaleAddOnLayoutWithoutCommunication(&ModRhs, &lambda, scale, m_InterMasterIndexLayout);
		}

	//	subtract solution on other processes from own value on gamma
		void compute_difference_on_gamma(vector_type& x)
		{
			VecSubtractOnLayout(&x,
								m_InterMasterIndexLayout,
								m_InterSlaveIndexLayout,
								&m_InterCommunicator);
			number scale = 1.0;
			if(pcl::GetProcRank() != 1)
				scale = -1.0;

			x *= scale;
		}

		void copy_dirichlet_values_and_zero(vector_type& ModRhs, const vector_type& r)
		{
			ModRhs.set(0.0);
			VecScaleAddOnLayoutWithoutCommunication(&ModRhs, &r, 1.0, m_InterSlaveIndexLayout);
			VecScaleAddOnLayoutWithoutCommunication(&ModRhs, &r, 1.0, m_InterMasterIndexLayout);
		}

		void set_dirichlet_rows_on_gamma(matrix_type& mat)
		{
			MatSetDirichletWithoutCommunication(&mat, m_InterSlaveIndexLayout);
			MatSetDirichletWithoutCommunication(&mat, m_InterMasterIndexLayout);
		}

	protected:
	//	set layouts to feti - partition local interchange
		void set_layouts_for_local_feti_partition_communication(vector_type& u)
		{
			u.set_slave_layout(m_LocalSlaveIndexLayout);
			u.set_master_layout(m_LocalMasterIndexLayout);
			u.set_communicator(m_LocalCommunicator);
			u.set_process_communicator(m_LocalProcessCommunicator);
			// todo: set vertical layouts iff gmg used
		}

	//	set layouts to inter-feti-partition interchange
		void set_layouts_for_inter_feti_partition_communication(vector_type& u)
		{
			u.set_slave_layout(m_InterSlaveIndexLayout);
			u.set_master_layout(m_InterMasterIndexLayout);
			u.set_communicator(m_InterCommunicator);
			u.set_process_communicator(m_InterProcessCommunicator);
			// todo: set vertical layouts iff gmg used
		}

	//	prepare local and inter feti-partition interfaces
		void prepare_layouts_and_communicators(vector_type& u)
		{
			m_InterSlaveIndexLayout = u.get_slave_layout();
			m_InterMasterIndexLayout = u.get_master_layout();
			m_InterCommunicator = u.get_communicator();
			m_InterProcessCommunicator = u.get_process_communicator();
			// todo: generalize to more than one process per feti-partition
		}

	protected:
	//	Feti-Partition-Local Interfaces and Communicators (i.e. non-gamma boundaries)
		IndexLayout m_LocalSlaveIndexLayout;
		IndexLayout m_LocalMasterIndexLayout;
		pcl::ProcessCommunicator m_LocalProcessCommunicator;
		pcl::ParallelCommunicator<IndexLayout> m_LocalCommunicator;

	//	Inter-Feti-Partition Interfaces and Communicators (i.e. for Gamma Boundary)
		IndexLayout m_InterMasterIndexLayout;
		IndexLayout m_InterSlaveIndexLayout;
		pcl::ProcessCommunicator m_InterProcessCommunicator;
		pcl::ParallelCommunicator<IndexLayout> m_InterCommunicator;

	protected:
	//	scaling of correction of flux
		number m_theta;

	// 	Operator that is inverted by this Inverse Operator
		IMatrixOperator<vector_type,vector_type,matrix_type>* m_A;

	// 	Parallel Matrix to invert
		matrix_type* m_pMatrix;

	//	Copy of matrix
		PureMatrixOperator<vector_type, vector_type, matrix_type> m_DirichletOperator;

	// 	Linear Solver to invert the local neumann problems
		ILinearOperatorInverse<vector_type,vector_type>* m_pNeumannSolver;

	// 	Linear Solver to invert the local dirichlet problems
		ILinearOperatorInverse<vector_type,vector_type>* m_pDirichletSolver;

	// 	Convergence Check
		IConvergenceCheck* m_pConvCheck;
};

#endif

} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__FETI__ */
