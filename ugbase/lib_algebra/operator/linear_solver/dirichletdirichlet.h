/*
 * dirichletdirichlet.h
 *
 *  Created on: 11.11.2010
 *      Author: iheppner, avogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__DIRICHLETDIRICHLET__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__DIRICHLETDIRICHLET__

namespace ug{

#ifdef UG_PARALLEL

#include <iostream>
#include <sstream>
#include <string>
#include "lib_algebra/operator/operator_inverse_interface.h"
#include "lib_algebra/parallelization/parallelization.h"
#include "lib_algebra/operator/debug_writer.h"

/// DirichletDirichletSolver implements (to be honest only) a Dirichlet-Dirichlet-solver.
/**
 * DirichletDirichletSolver implements a Dirichlet-Dirichlet-solver for a partitioning
 * into two nonoverlapping subdomains (e.g. "FETI subdomains"). See e.g.
 * "Domain Decomposition Methods -- Algorithms and Theory",
 * A. Toselli, O. Widlund, Springer 2004, sec. 1.3.5, p. 12ff.
 *
 * REMARKS:
 * We have to consider two types of interfaces/layouts, resp. communication
 * over these interfaces here:
 *
 *    + between "FETI subdomains", in general sub-partitioned into several
 *      "processor partitions" (i.e.: number of procs per FETI subdomain > 1),
 *      over the "FETI-Gamma"
 *      ==> "inter" (subdomain) communication (formerly also called "global" communication),
 *
 *      and
 *    + between "processor partitions" inside an FETI subdomain,
 *      over the non "FETI-Gamma"
 *      ==> "intra" subdomain communication (formerly also called "local", "inner" communication).
 *
 * \TAlgebra	type of algebra (template parameter), e.g. CPUAlgebra.
 */
template <typename TAlgebra>
class DirichletDirichletSolver : public IMatrixOperatorInverse<	typename TAlgebra::vector_type,
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
		DirichletDirichletSolver() :
			m_theta(1.0), m_pOperator(NULL), m_pNeumannSolver(NULL),
			m_pDirichletSolver(NULL), m_pConvCheck(NULL),
			m_pDebugWriter(NULL)
		{}

	///	name of solver
		virtual const char* name() const {return "Dirichlet-Dirichlet Solver";}

	///	sets a convergence check
		void set_convergence_check(IConvergenceCheck& convCheck)
		{
			m_pConvCheck = &convCheck;
			m_pConvCheck->set_offset(3);
		}

	/// returns the convergence check
		IConvergenceCheck* get_convergence_check() {return m_pConvCheck;}

	///	sets a sequential Neumann solver
		void set_neumann_solver(ILinearOperatorInverse<vector_type, vector_type>& neumannSolver)
		{
			m_pNeumannSolver = &neumannSolver;
		}

	///	sets a sequential Dirichlet solver
		void set_dirichlet_solver(ILinearOperatorInverse<vector_type, vector_type>& dirichletSolver)
		{
			m_pDirichletSolver = &dirichletSolver;
		}

	///	sets damping factor
		void set_theta(number theta)
		{
			m_theta = theta;
		}

	//	set debug output
		void set_debug(IDebugWriter<algebra_type>* debugWriter)
		{
			m_pDebugWriter = debugWriter;
		}

	///	initializes the solver for operator A
		virtual bool init(MatrixOperator<vector_type, vector_type, matrix_type>& A)
		{
		//	save current operator
			m_pOperator = &A;

		//	save matrix to invert
			m_pMatrix = &m_pOperator->get_matrix();

		//	Copy Matrix for Dirichlet Problem
			matrix_type& dirMat = m_DirichletOperator.get_matrix();
			dirMat = *m_pMatrix;

		//	Set Dirichlet values
			set_dirichlet_rows_on_gamma(dirMat);

		//	init sequential solver for Neumann problem
			if(m_pNeumannSolver != NULL)
				if(!m_pNeumannSolver->init(*m_pOperator))
				{
					UG_LOG("ERROR in 'DirichletDirichletSolver::init': Cannot init "
							"Sequential Neumann Solver for Operator A.\n");return false;
				}

		//	init sequential solver for Dirichlet problem
			if(m_pDirichletSolver != NULL)
				if(!m_pDirichletSolver->init(m_DirichletOperator))
				{
					UG_LOG("ERROR in 'DirichletDirichletSolver::init': Cannot init "
							"Sequential Dirichlet Solver for Operator A.\n");return false;
				}

		//	check that solvers are different
			if(m_pNeumannSolver == m_pDirichletSolver)
			{
				UG_LOG("ERROR in 'DirichletDirichletSolver:prepare': Solver for dirichlet"
						" and neumann problem must be different instances.\n");
				return false;
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
		virtual bool apply_return_defect(vector_type& x, vector_type& b)
		{
			UG_LOG("INFO: parameter 'theta' used in 'DirichletDirichletSolver::apply_return_defect' is '" << m_theta << "' (TMP)\n");

		//	check that matrix has been set
			if(m_pOperator == NULL)
			{
				UG_LOG("ERROR: In 'DirichletDirichletSolver::apply_return_defect': Matrix A not set.\n");
				return false;
			}

		//	check that sequential solver has been set
			if(m_pNeumannSolver == NULL)
			{
				UG_LOG("ERROR: In 'DirichletDirichletSolver::apply_return_defect': No sequential Neumann Solver set.\n");
				return false;
			}
			if(m_pDirichletSolver == NULL)
			{
				UG_LOG("ERROR: In 'DirichletDirichletSolver::apply_return_defect': No sequential Dirichlet Solver set.\n");
				return false;
			}

		//	check that convergence check is set
			if(m_pConvCheck == NULL)
			{
				UG_LOG("ERROR: In 'DirichletDirichletSolver::apply_return_defect':"
						" Convergence check not set.\n");
				return false;
			}

		//	Check parallel storage type of vectors
			if(!b.has_storage_type(PST_ADDITIVE) || !x.has_storage_type(PST_CONSISTENT))
			{
				UG_LOG("ERROR: In 'DirichletDirichletSolver::apply_return_defect': "
						"Inadequate storage format of Vectors.\n");
				return false;
			}

		// 	todo: 	it would be sufficient to only copy the pattern (and parallel constructor)
		//			without initializing the values

		//	Normal Flux at FETI-Gamma (i.e. inner boundary of feti method)
			vector_type lambda; lambda.create(x.size()); lambda = x;

		//	Flux Correction at Gamma - Boundary
			vector_type eta; eta.create(x.size()); eta = x;

		//	Modified right - hand side for Neumann problem
			vector_type ModRhs;         ModRhs.create(b.size()); ModRhs     = b;
			vector_type ModRhsCopy; ModRhsCopy.create(b.size()); ModRhsCopy = b;

		//  Prepare convergence check
			prepare_conv_check();

		//	flag iff first iterate
			bool first = true;

		//	set initial lambda (flux over \Gamma) to zero vector
			lambda.set(0.0);

		// 	Iteration loop
			m_iterCnt = 0;
			while(true)
			{
				m_iterCnt++;

			//	Copy right-hand side
				ModRhs = b; // we always add the new flux iterate to the *original* rhs!

			//	set inter subdomain communication (added 07122010ih) - always before "FETI operations"
			// (hier und nicht erst am Schluss der Schleife, evt. wichtig fuer andere Startwerte von lambda!)
				//x.use_layout(1);
				//ModRhs.use_layout(1);


			//	write debug
				write_debug(lambda, "FetiLambdaOnGamma");

			//	Compute new rhs for Neumann problem using lambda on ("FETI-") Gamma
				add_flux_to_rhs(ModRhs, lambda);

			//	todo: Here and in the following, adapt to Bnd-Layout
			//	set intra subdomain communication - always before solution step
				//x.use_layout(0);
				//ModRhs.use_layout(0);

			//	Solve Neumann problem on FETI subdomain
				if(!m_pNeumannSolver->apply_return_defect(x, ModRhs))
				{
					UG_LOG_ALL_PROCS("ERROR in 'DirichletDirichletSolver::apply_return_defect': "
									 "Could not solve Neumann problem on proc " << pcl::GetProcRank() << ".\n");
					return false;
				}
				// todo: Check that all processes solved the problem

			//	write debug
				write_debug(x, "FetiNeumann");

			//	set inter subdomain communication
				//x.use_layout(1);
				//ModRhs.use_layout(1);

			//	Set Dirichlet values for Rhs, zero else
				copy_dirichlet_values_and_zero(ModRhs, x);

			//	Compute difference of solution on Gamma
				compute_difference_on_gamma(ModRhs);

			//	write debug
				write_debug(ModRhs, "FetiDiffSolOnGamma");

				ModRhsCopy = ModRhs; /* 'ModRhsCopy' gets also layouts of 'ModRhs's (22122010ih) */

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

			//	set intra subdomain communication
				//x.use_layout(0);
				//ModRhs.use_layout(0);

			//	Solve Dirichlet problem on FETI subdomain
				if(!m_pDirichletSolver->apply_return_defect(x, ModRhs))
				{
					UG_LOG_ALL_PROCS("ERROR in 'DirichletDirichletSolver::apply_return_defect': "
									 "Could not solve Dirichlet problem on proc " << pcl::GetProcRank() << ".\n");
					return false;
				}
				// todo: Check that all processes solved the problem

			//	write debug
				write_debug(x, "FetiDirichlet");

			//	set inter subdomain communication
				//x.use_layout(1);
				//eta.use_layout(1);

			//	Compute update for lambda: eta = A*x (to be more specific: \eta_i^{n+1} = A_{\Gamma I}^{(i)} w_i^{n+1} + A_{\Gamma \Gamma}^{(i)} r_{\Gamma})
				m_pMatrix->apply(eta, x);

			//	sum up over all processes
				eta.set_storage_type(PST_ADDITIVE);
				eta.change_storage_type(PST_CONSISTENT);

			// 	Update lambda
				VecScaleAdd(lambda, 1.0, lambda, -m_theta, eta);

			} /* end of iteration loop */

		//	Post Output
			if(!m_pConvCheck->post())
			{
				UG_LOG("ERROR in 'DirichletDirichletSolver::apply_return_defect': "
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
		virtual ~DirichletDirichletSolver() {};

	protected:
	//	Prepare the convergence check
		void prepare_conv_check()
		{
			m_pConvCheck->set_name(name());
			m_pConvCheck->set_symbol('%');
			m_pConvCheck->set_name(name());

			if(m_pNeumannSolver != NULL && m_pDirichletSolver != NULL)
			{
				std::stringstream ss; ss <<  " (Seq. Solver: " << m_pNeumannSolver->name() << ","
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
			/* mrupp: removed warning unused variable scale
			number scale = 1.0;
			if(pcl::GetProcRank() != 0) // TODO: generalize to more than one process per FETI subdomain 
					scale = -1.0;
			*/

		// TODO: Implement Bnd-Layout

		//	VecScaleAddOnLayout(&ModRhs, &lambda, scale, ModRhs.slave_layout(1));
		//	VecScaleAddOnLayout(&ModRhs, &lambda, scale, ModRhs.master_layout(1));
		}

	//	subtract solution on other processes from own value on gamma
		void compute_difference_on_gamma(vector_type& x)
		{
		// TODO: Implement Bnd-Layout
		/* VecSubtractOnLayout(&x,
								x.master_layout(1),
								x.slave_layout(1),
								&x.communicator(1));*/
			number scale = 1.0;
			if(pcl::GetProcRank() != 0) // TODO: generalize to more than one process per FETI subdomain
				scale = -1.0;

			x *= scale;
		}

		void copy_dirichlet_values_and_zero(vector_type& ModRhs, const vector_type& r)
		{
			ModRhs.set(0.0);
		//	todo: implement Bnd-Layout
		//	VecScaleAddOnLayout(&ModRhs, &r, 1.0, ModRhs.slave_layout(1));
		//	VecScaleAddOnLayout(&ModRhs, &r, 1.0, ModRhs.master_layout(1));
		}

		void set_dirichlet_rows_on_gamma(matrix_type& mat)
		{
		//	todo: implement Bnd-Layout
		//	MatSetDirichletOnLayout(&mat, mat.slave_layout(1));
		//	MatSetDirichletOnLayout(&mat, mat.master_layout(1));
		}

	protected:
		void write_debug(const vector_type& vec, const char* filename)
		{
		//	add iter count to name
			std::string name(filename);
			char ext[20]; sprintf(ext, "_iter%03d", m_iterCnt);
			name.append(ext);

		//	if no debug writer set, we're done
			if(m_pDebugWriter == NULL) return;

		//	write
			m_pDebugWriter->write_vector(vec, name.c_str());
		}

		int m_iterCnt;

	protected:
	//	scaling of correction of flux
		number m_theta;

	// 	Operator that is inverted by this Inverse Operator
		MatrixOperator<vector_type,vector_type,matrix_type>* m_pOperator;

	// 	Parallel Matrix to invert
		matrix_type* m_pMatrix;

	//	Copy of matrix
		MatrixOperator<vector_type, vector_type, matrix_type> m_DirichletOperator;

	// 	Linear Solver to invert the local Neumann problems
		ILinearOperatorInverse<vector_type,vector_type>* m_pNeumannSolver;

	// 	Linear Solver to invert the local Dirichlet problems
		ILinearOperatorInverse<vector_type,vector_type>* m_pDirichletSolver;

	// 	Convergence Check
		IConvergenceCheck* m_pConvCheck;

	//	Debug Writer
		IDebugWriter<algebra_type>* m_pDebugWriter;
};

#endif

} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__DIRICHLETDIRICHLET__ */
