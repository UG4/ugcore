// 01.02.2011 (m,d,y)

//	THIS FILE IS ONLY TEMPORARY!
//	I (Sebastian) created this file to reduce compile-time during
//	the creation of the FETI solver. When everything works as
//	expected one should move the code in this file to
//	feti_impl.h, so that it will work with all the different algebras.
//
//	In the moment template instantiations are invoked at the end of this file.

//	do not include feti.h itself or your will be sent to compile-error-hell
//#include "feti.h"
#include "lib_algebra/lib_algebra.h"
#include <cmath>

namespace ug{

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	LocalSchurComplement implementation
template <typename TAlgebra>
LocalSchurComplement<TAlgebra>::
LocalSchurComplement() :
	m_pMatrix(NULL),
	m_pMasterPrimalLayout(NULL),
	m_pSlavePrimalLayout(NULL),
	m_pSlaveDualLayout(NULL),
	m_pMasterDualLayout(NULL),
	m_pSlaveDualNbrLayout(NULL),
	m_pMasterDualNbrLayout(NULL),
	m_pDirichletMatrix(NULL),
	m_pDirichletSolver(NULL),
	m_pDebugWriter(NULL)
{
}

template <typename TAlgebra>
bool LocalSchurComplement<TAlgebra>::
init()
{
//	check that operator has been set
	if(m_pOperator == NULL)
	{
		UG_LOG("ERROR in 'LocalSchurComplement::init': No Operator A"
				" set.\n");
		return false;
	}

//	check that Primal layouts have been set
	if(m_pSlavePrimalLayout == NULL || m_pMasterPrimalLayout == NULL)
	{
		UG_LOG("ERROR in 'LocalSchurComplement::init': Master or Slave"
				" layout for primal unknowns not set.\n");
		return false;
	}

//	check that Primal layouts have been set
	if(m_pSlaveDualLayout == NULL || m_pMasterDualLayout == NULL)
	{
		UG_LOG("ERROR in 'LocalSchurComplement::init': Master or Slave"
				" layout for dual unknowns not set.\n");
		return false;
	}

//	save matrix from which we build the Schur complement
	m_pMatrix = &m_pOperator->get_matrix();

//	get matrix from dirichlet operator
	m_pDirichletMatrix = &m_DirichletOperator.get_matrix();

//	Copy Matrix for Dirichlet Problem
	*m_pDirichletMatrix = *m_pMatrix;

//	Set Dirichlet values on Pi
	MatSetDirichletOnLayout(m_pDirichletMatrix, *m_pSlavePrimalLayout);
	MatSetDirichletOnLayout(m_pDirichletMatrix, *m_pMasterPrimalLayout);

//	Set Dirichlet values on Delta
	MatSetDirichletOnLayout(m_pDirichletMatrix, *m_pSlaveDualLayout);
	MatSetDirichletOnLayout(m_pDirichletMatrix, *m_pMasterDualLayout);

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
		/*m_pDebugWriter->write_matrix(m_DirichletOperator.get_matrix(),
									 "FetiDirichletMatrix");
		m_pDebugWriter->write_matrix(m_pOperator->get_matrix(),
									 "FetiOriginalMatrix");*/
	}

//	we're done
	return true;
} /* end 'LocalSchurComplement::init()' */

template <typename TAlgebra>
bool LocalSchurComplement<TAlgebra>::
apply(vector_type& f, const vector_type& u)
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

//	Check parallel storage type of matrix
	if(!m_pDirichletMatrix->has_storage_type(PST_ADDITIVE))
	{
		UG_LOG("ERROR: In 'LocalSchurComplement::apply': "
						"Inadequate storage format of matrix.\n");
		return false;
	}

//	Check parallel storage type of vectors
	if (!u.has_storage_type(PST_CONSISTENT))
	{
		UG_LOG("ERROR: In 'LocalSchurComplement::apply': "
						"Inadequate storage format of Vector 'u' (should be consistent).\n");
		return false;
	}
	if(!f.has_storage_type(PST_ADDITIVE))
	{
		UG_LOG("ERROR: In 'LocalSchurComplement::apply': "
						"Inadequate storage format of Vector 'f' (should be additive).\n");
		return false;
	}

//	Help vector
//	\todo: it would be sufficient to copy only the layouts without copying the values
	vector_type uTmp; uTmp.create(u.size()); uTmp = u;

//	1. Set values to zero on I and Pi
	// (a) Reset all values
	uTmp.set(0.0);

	// (b) Copy values on \Delta
	VecScaleAppendOnLayout(&uTmp, &u, 1.0, *m_pSlaveDualLayout);
	VecScaleAppendOnLayout(&uTmp, &u, 1.0, *m_pMasterDualLayout);

//	2. Compute rhs f_{I} = A_{I \Delta} u_{\Delta}
	if(!m_DirichletOperator.apply(f, uTmp))
	{
		UG_LOG_ALL_PROCS("ERROR in 'LocalSchurComplement::apply': "
						 "Could not compute Rhs for Dirichlet problem on "
						 "Proc " << pcl::GetProcRank() << ".\n");
		return false;
	}
	// set values to zero on \Delta
	VecSetOnLayout(&f, 0.0, *m_pSlaveDualLayout);
	VecSetOnLayout(&f, 0.0, *m_pMasterDualLayout);

//	3. Invert on inner unknowns u_{I} = A_{II}^{-1} f_{I}
	// (a) use the inner-FETI-block layouts
	//uTmp.use_layout(0);
	//f.use_layout(0);
	//m_pDirichletMatrix->use_layout(0);

	// (b) invoke Dirichlet solver
	if(!m_pDirichletSolver->apply_return_defect(uTmp, f))
	{
		UG_LOG_ALL_PROCS("ERROR in 'LocalSchurComplement::apply': "
						 "Could not solve Dirichlet problem on Proc "
							<< pcl::GetProcRank() << ".\n");
		return false;
	}

//	4. Compute result vector
	// (a) Scale u_{I} by -1
	uTmp *= -1.0;

	// (b) Add u_{\Delta} on \Delta
	VecScaleAppendOnLayout(&uTmp, &u, 1.0, *m_pSlaveDualLayout);
	VecScaleAppendOnLayout(&uTmp, &u, 1.0, *m_pMasterDualLayout);

	// (c) Multiply with full matrix
	if(!m_pOperator->apply(f, uTmp))
	{
		UG_LOG_ALL_PROCS("ERROR in 'LocalSchurComplement::apply': "
						 "Could not apply full matrix on "
						 "Proc " << pcl::GetProcRank() << ".\n");
		return false;
	}

//	5. Reset all values for I, \Pi
	VecSetExcludingLayout(&f, 0.0, *m_pSlaveDualLayout);

	VecSetOnLayout(&f, 0.0, *m_pSlavePrimalLayout);
	VecSetOnLayout(&f, 0.0, *m_pMasterPrimalLayout);

//	we're done
	return true;
} /* end 'LocalSchurComplement::apply()' */

template <typename TAlgebra>
bool LocalSchurComplement<TAlgebra>::
apply_sub(vector_type& f, const vector_type& u)
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

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	SchurComplementInverse implementation
template <typename TAlgebra>
SchurComplementInverse<TAlgebra>::
SchurComplementInverse() :
	m_A(NULL),
	m_pMatrix(NULL),
	m_pMasterPrimalLayout(NULL),
	m_pSlavePrimalLayout(NULL),
	m_pSlaveDualLayout(NULL),
	m_pMasterDualLayout(NULL),
	m_pSlaveDualNbrLayout(NULL),
	m_pMasterDualNbrLayout(NULL),
	m_pNeumannMatrix(NULL),
	m_pNeumannSolver(NULL),
	m_primalRootProc(-1),
	m_pOneProcSchurCompMatrix(NULL),
	m_pConvCheck(NULL),
	m_pDebugWriter(NULL)
{
}

template <typename TAlgebra>
bool SchurComplementInverse<TAlgebra>::
init(ILinearOperator<vector_type, vector_type>& L)
{
//	success flag
	bool bSuccess = true;

//	remember operator
	m_A = dynamic_cast<IMatrixOperator<vector_type, vector_type, matrix_type>*>(&L);

//	check, that operator is correct
	if(m_A == NULL)
	{
		UG_LOG_ALL_PROCS("ERROR in 'SchurComplementInverse::init':"
				" Wrong type of operator passed for init.\n");
		bSuccess = false;
	}

//	check that Pi layouts have been set
	if(m_pSlavePrimalLayout == NULL || m_pMasterPrimalLayout == NULL)
	{
		UG_LOG_ALL_PROCS("ERROR in 'SchurComplementInverse::init':"
				" Master or Slave layout for cross points not set "
				"on Proc " << pcl::GetProcRank() << ".\n");
		bSuccess = false;
	}

//	Check all procs
	if(!pcl::AllProcsTrue(bSuccess))
	{
		UG_LOG_ALL_PROCS("ERROR in 'SchurComplementInverse::init':"
				" Some proc could not init Schur Complement inverse.\n");
		return false;
	}

//	save matrix from which we build the Schur complement
	m_pMatrix = &m_A->get_matrix();

//	get matrix from Neumann operator
	m_pNeumannMatrix = &m_NeumannOperator.get_matrix();

//	Copy Matrix for Neumann Problem
	*m_pNeumannMatrix = *m_pMatrix;

//	Set Dirichlet values on Pi
	MatSetDirichletOnLayout(m_pNeumannMatrix, *m_pSlavePrimalLayout);
	MatSetDirichletOnLayout(m_pNeumannMatrix, *m_pMasterPrimalLayout);

//	init sequential solver for Dirichlet problem
	if(m_pNeumannSolver != NULL)
		if(!m_pNeumannSolver->init(m_NeumannOperator))
		{
			UG_LOG("ERROR in 'SchurComplementInverse::init': Cannot init "
					"Sequential Neumann Solver for Operator A.\n");
			return false;
		}

//	Debug output of matrices
	if(m_pDebugWriter != NULL)
	{
		m_pDebugWriter->write_matrix(m_NeumannOperator.get_matrix(),
									 "FetiNeumannMatrix");
	}

//	Choose root process, where Schur complement w.r.t. Primal unknowns
//	is gathered.
	m_primalRootProc = pcl::GetOutputProcRank();

//	vector to store newly created root ids
	std::vector<int> newMasterIDsOut;

//	Build layouts such that all processes can communicated their unknowns
//	to the primal Root Process
	int newVecSize = BuildOneToManyLayout(m_masterAllToOneLayout,
						 m_slaveAllToOneLayout, m_primalRootProc,
						 *m_pMasterPrimalLayout, *m_pSlavePrimalLayout,
						 m_allToOneProcessComm, &newMasterIDsOut);

//	We have to gather the quantities of primal nodes on each process
//	of the feti-block in one array.
	int numLocalPrimals = m_slaveAllToOneLayout.num_interface_elements();
	UG_LOG("num local primals: " << numLocalPrimals << endl);
	m_primalQuantities.resize(m_localFetiBlockComm.size());
	m_localFetiBlockComm.allgather(&numLocalPrimals, 1, PCL_DT_INT,
							&m_primalQuantities.front(), 1, PCL_DT_INT);

//	log num primal quantities
	UG_LOG("proc ids: ");
	for(size_t i = 0; i < m_primalQuantities.size(); ++i){
		UG_LOG(m_localFetiBlockComm.get_proc_id(i) << " ");
	}
	UG_LOG(endl);
	UG_LOG("primal quantities: ");
	for(size_t i = 0; i < m_primalQuantities.size(); ++i){
		UG_LOG(m_primalQuantities[i] << " ");
	}
	UG_LOG(endl);

//	build matrix on primalRoot
	if(pcl::GetProcRank() == m_primalRootProc)
	{
	//	get matrix
		m_pOneProcSchurCompMatrix = &m_OneProcSchurCompOp.get_matrix();

	//	create matrix of correct size
		m_pOneProcSchurCompMatrix->create(newVecSize, newVecSize);

		std::cout << "On PrimalRoot: Creating proc local Schur Complement"
					" of size " << newVecSize <<"x"<<newVecSize << std::endl;
	}

//	we're done
	return true;
} /* end 'SchurComplementInverse::init()' */



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	FETISolver implementation
template <typename TAlgebra>
FETISolver<TAlgebra>::
FETISolver() :
	m_A(NULL),
	m_pMatrix(NULL),
	m_pDirichletSolver(NULL),
	m_pNeumannSolver(NULL),
	m_pConvCheck(NULL),
	m_pDebugWriter(NULL)
{

}

template <typename TAlgebra>
bool FETISolver<TAlgebra>::
init(IMatrixOperator<vector_type, vector_type, matrix_type>& A)
{
//	bool flag
	bool bSuccess = true;

//	remember A
	m_A = &A;

//	get matrix
	m_pMatrix = &m_A->get_matrix();

//	create FETI Layouts:
// 	\todo: @Ingo, Doc Me please
	BuildDomainDecompositionLayouts(m_masterDualLayout, m_slaveDualLayout,
					m_masterInnerLayout, m_slaveInnerLayout, m_masterDualNbrLayout,
					m_slaveDualNbrLayout, m_masterPrimalLayout, m_slavePrimalLayout,
					m_pMatrix->get_master_layout(), m_pMatrix->get_slave_layout(),
					(int)(m_pMatrix->num_rows() - 1), *m_pDDInfo);

//	create local feti block communicator
	int localSubdomID = m_pDDInfo->map_proc_id_to_subdomain_id(pcl::GetProcRank());
	pcl::ProcessCommunicator worldComm;
	for(int i = 0; i < m_pDDInfo->get_num_subdomains(); ++i){
		if(localSubdomID == i)
			m_localFetiBlockComm = worldComm.create_sub_communicator(true);
		else
			worldComm.create_sub_communicator(false);
	}

	pcl::ParallelCommunicator<IndexLayout> comTmp;
	PrintLayout(comTmp, m_masterInnerLayout, m_slaveInnerLayout);
/*
//	write layouts
	pcl::SynchronizeProcesses();
	UG_LOG("------------- DUAL MASTER ------------\n")
	LogIndexLayoutOnAllProcs(m_masterDualLayout, 1);
	pcl::SynchronizeProcesses();
	UG_LOG("------------- DUAL SLAVE  ------------\n")
	LogIndexLayoutOnAllProcs(m_slaveDualLayout, 1);

	pcl::SynchronizeProcesses();
	UG_LOG("------------- PRIMAL MASTER ------------\n")
	LogIndexLayoutOnAllProcs(m_masterPrimalLayout, 1);
	pcl::SynchronizeProcesses();
	UG_LOG("------------- PRIMAL SLAVE  ------------\n")
	LogIndexLayoutOnAllProcs(m_slavePrimalLayout, 1);
*/
//  ----- INIT DIRICHLET SOLVER  ----- //

//	check that dirichlet solver has been set
	if(m_pDirichletSolver == NULL)
	{
		UG_LOG("ERROR in FETISolver::init: No dirichlet solver set "
				" for inversion of A_{II} in Local Schur complement.\n");
		return false;
	}

//	set layouts
	m_LocalSchurComplement.set_primal_layouts(m_slavePrimalLayout, m_masterPrimalLayout);
	m_LocalSchurComplement.set_dual_layouts(m_slaveDualLayout, m_masterDualLayout);

//	set dirichlet solver for local schur complement
	m_LocalSchurComplement.set_dirichlet_solver(*m_pDirichletSolver);

//	set operator in local schur complement
	m_LocalSchurComplement.set_matrix(*m_A);

//	init local Schur complement
	if(m_LocalSchurComplement.init() != true)
	{
		UG_LOG("ERROR in FETISolver::init: Can not init local Schur "
				"complement.\n");
		bSuccess = false;
	}

//	check all procs
	if(!pcl::AllProcsTrue(bSuccess))
	{
		UG_LOG("ERROR in FETISolver::init: Some process could not init"
				" local Schur complement.\n");
		return false;
	}

//  ----- INIT NEUMANN SOLVER  ----- //

//	check that neumann solver has been set
	if(m_pNeumannSolver == NULL)
	{
		UG_LOG("ERROR in FETISolver::init: No neumann solver set "
				" for inversion of A_{I,Delta}{I,Delta} in Local Schur complement.\n");
		return false;
	}

//	set the local feti block communicator
	m_SchurComplementInverse.set_local_feti_block_communicator(m_localFetiBlockComm);

//	set layouts
	m_SchurComplementInverse.set_primal_layouts(m_slavePrimalLayout, m_masterPrimalLayout);
	m_SchurComplementInverse.set_dual_layouts(m_slaveDualLayout, m_masterDualLayout);

//	set neumann solver in SchurComplementInverse
	m_SchurComplementInverse.set_neumann_solver(*m_pNeumannSolver);

//	init Schur complement inverse
	if(m_SchurComplementInverse.init(*m_A) != true)
	{
		UG_LOG("ERROR in FETISolver::init: Can not init Schur "
				"complement inverse.\n");
		bSuccess = false;
	}

//	check all procs
	if(!pcl::AllProcsTrue(bSuccess))
	{
		UG_LOG("ERROR in FETISolver::init: Some process could not init"
				" Schur complement inverse.\n");
		return false;
	}

//  ----- INIT DUAL DIRICHLET SOLVER  ----- //

//	check that dual dirichlet solver has been set
	if(m_pDualDirichletSolver == NULL)
	{
		UG_LOG("ERROR in FETISolver::init: No dual dirichlet solver set "
				" for computation of tilde{f}_Delta.\n");
		return false;
	}

//	get matrix from dual dirichlet operator
	m_pDualDirichletMatrix = &m_DualDirichletOperator.get_matrix();

//	Copy Matrix
	*m_pDualDirichletMatrix = *m_pMatrix;

//	Set Dirichlet values on Dual
	MatSetDirichletOnLayout(m_pDualDirichletMatrix, m_slaveDualLayout);
	MatSetDirichletOnLayout(m_pDualDirichletMatrix, m_masterDualLayout);

//	make diagonal of matrix consistent on Pi
	MatAdditiveToConsistentOnDiag<algebra_type>(m_pDualDirichletMatrix, m_masterDualLayout,
	                              m_slaveDualLayout);

//	init sequential solver for Dirichlet problem
	if(m_pDualDirichletSolver != NULL)
		if(!m_pDualDirichletSolver->init(m_DualDirichletOperator))
		{
			UG_LOG("ERROR in 'SchurComplementInverse::init': Cannot init "
					"Dual Dirichlet Solver for Operator A.\n");
			return false;
		}

//	we're done
	return true;
} /* end 'FETISolver::init()' */

template <typename TAlgebra>
bool FETISolver<TAlgebra>::
apply_return_defect(vector_type& u, vector_type& f)
{
//	This function is used to solve the system Au=f. While the matrix A has
//	already been passed in additive storage in the function init(), here are
//	passed the (unknown, to be computed) solution u, that must be returned in
//	consistent storage, and the right-hand side, that is given in additive
//	storage. Note, that the right-hand side is overwritten, such that it can
//	be used for the computation of the defect (Once the defect is computed the
//	right-hand side is not needed anymore).

//	check, that convergence check has been set
	if(m_pConvCheck == NULL)
	{
		UG_LOG("ERROR: In 'FETISolver::apply': Convergence check not set.\n");
		return false;
	}

//	check storage type
	if(!f.has_storage_type(PST_ADDITIVE) || !u.has_storage_type(PST_CONSISTENT))
	{
		UG_LOG("ERROR: In 'FETISolver::apply':Inadequate storage format of "
				"Vectors 'u' and 'f'.\n");
		return false;
	}

//	Construct some vectors, that are all needed on the Dual unknowns.
//	These vectors are used exclusively on the Dual unknowns, but for facility we
//	use storage for a whole vector and only use it Dual entries.

//	lagrange multiplier
	vector_type lambda; lambda.create(u.size());

//	residuum
	vector_type r; r.create(u.size());

//	search direction
	vector_type p; p.create(u.size());

//	preconditioned residuum
	vector_type z; z.create(u.size());

//	help vector to compute t = F*p
	vector_type t; t.create(u.size());


//	reset start value of lagrange multiplier
	lambda.set(0.0);


//	The idea of the Feti solver is to perform a cg-method on the Dual unknowns.
//	Therefore, the right-hand side has to be adapted, such that f_Dual fits to
//	the Schur complement system on the Dual unknowns. i.e.
//	\tilde{f}_{\Delta} = f_{\Delta}
//						 - A_{\Delta \{I \Pi\}} (A_{\{I \Pi\} \{I \Pi\}})^{-1} f_{\{I \Pi\}}
//
//	Please note: To compute this, we have to make the matrix consistent in the
//				 Primal unknowns. Then, the rhs can be computed on each process
//				 individually.

//	make f consistent on Pi
	AdditiveToConsistent(&f, m_masterPrimalLayout, m_slavePrimalLayout);

//	set inner layouts
	t.set_slave_layout(m_slaveInnerLayout);
	t.set_master_layout(m_masterInnerLayout);
	t.set_process_communicator(m_localFetiBlockComm);

	f.set_slave_layout(m_slaveInnerLayout);
	f.set_master_layout(m_masterInnerLayout);
	f.set_process_communicator(m_localFetiBlockComm);

//	compute \tilde{f}_{\Delta}
	if(!compute_tilde_f(t, f)) return false;

// 	Build start residuum:  r = r0 := d - F*lambda.
//	This is done in three steps.

//	(a) Compute d:= B S_{\Delta \Delta}^{-1} \tilde{f}_{\Delta}
//		and set r0 = d
	if(!compute_d(r, t)) return false;

// 	(b) Build t = F*lambda (t is additive afterwards)
	if(!apply_F(t, lambda))
	{
		UG_LOG("ERROR in 'FETISolver::apply': Unable "
			   "to build t = F*p. Aborting.\n"); return false;
	}

// (c) Subtract values on \Delta, r0 = r0 - t
	VecScaleAppendOnDual(r, t, -1.0);

// 	Precondition the start defect: apply z = M^-1 * r
	if (!apply_M_inverse_with_identity_scaling(z, r))
	{
		UG_LOG("ERROR in 'FETISolver::apply': "
			   "Cannot apply preconditioner. Aborting.\n");
		return false;
	}

//	prepare appearance of conv check
	prepare_conv_check();

//	compute and set start defect
	m_pConvCheck->start_defect(VecNormOnDual(r));

//	start values
	number rho, rho_new, beta, alpha, alpha_denominator;
	rho = rho_new = beta = alpha = alpha_denominator = 0.0;

// 	start search direction
	p = z;

// 	start rho
	rho = VecProdOnDual(z, r);

//	reset iteration count
	m_iterCnt = 0;

// 	Iteration loop
	while(!m_pConvCheck->iteration_ended())
	{
	//	increase iteration count
		m_iterCnt++;

	// 	Build t = F*p (t is additive afterwards)
		if(!apply_F(t, p))
		{
			UG_LOG("ERROR in 'FETISolver::apply': Unable "
						"to build t = F*p. Aborting.\n"); return false;
		}

	// 	Compute alpha
		alpha_denominator = VecProdOnDual(t, p);

		alpha = rho/alpha_denominator;

	// 	Update lambda := lambda + alpha*p
		VecScaleAddOnDual(lambda, 1.0, lambda, alpha, p);

	// 	Update r := r - alpha*t
		VecScaleAddOnDual(r, 1.0, r, -alpha, t);

	// 	Compute new defect
		m_pConvCheck->update_defect(VecNormOnDual(r));

	// 	Preconditioning: apply z = M^-1 * r
		if (!apply_M_inverse_with_identity_scaling(z, r))
		{
			UG_LOG("ERROR in 'FETISolver::apply': "
				   "Cannot apply preconditioner. Aborting.\n"); return false;
		}

	// 	new rho
		rho_new = VecProdOnDual(z, r);

	// 	new beta
		beta = rho_new/rho;

	// 	new direction p:= beta*p + z
		VecScaleAddOnDual(p, beta, p, 1.0, z);

	// 	update rho
		rho = rho_new;
	} /* end iteration loop */

	return m_pConvCheck->post();
} /* end 'FETISolver::apply_return_defect()' */

template <typename TAlgebra>
number FETISolver<TAlgebra>::
VecNormOnDual(vector_type& vec)
{
//	forward to VecProc
	return sqrt(VecProdOnDual(vec, vec));
}

template <typename TAlgebra>
void FETISolver<TAlgebra>::
VecScaleAddOnDual(vector_type& vecDest,
                  number alpha1, const vector_type& vecSrc1,
                  number alpha2, const vector_type& vecSrc2)
{
	VecScaleAddOnLayout(&vecDest, alpha1, &vecSrc1, alpha2, &vecSrc2, m_slaveDualLayout);
	VecScaleAddOnLayout(&vecDest, alpha1, &vecSrc1, alpha2, &vecSrc2, m_masterDualLayout);
}

template <typename TAlgebra>
void FETISolver<TAlgebra>::
VecScaleAppendOnDual(vector_type& vecInOut,
                     const vector_type& vecSrc1, number alpha1)
{
	VecScaleAppendOnLayout(&vecInOut, &vecSrc1, alpha1, m_slaveDualLayout);
	VecScaleAppendOnLayout(&vecInOut, &vecSrc1, alpha1, m_masterDualLayout);
}

template <typename TAlgebra>
void FETISolver<TAlgebra>::
VecSetOnDual(vector_type& vecSrc, number alpha)
{
	VecSetOnLayout(&vecSrc, alpha, m_slaveDualLayout);
	VecSetOnLayout(&vecSrc, alpha, m_masterDualLayout);
}

template <typename TAlgebra>
number FETISolver<TAlgebra>::
VecProdOnDual(const vector_type& vecSrc1, const vector_type& vecSrc2)
{
//	reset result
	number prod = 0.0, prodTmp = 0.0;

//	add prod on master
	VecProdOnLayout(prodTmp, &vecSrc1, &vecSrc2, m_masterDualLayout);
	prod += prodTmp;

//	add prod on slave
	VecProdOnLayout(prodTmp, &vecSrc1, &vecSrc2, m_slaveDualLayout);
	prod += prodTmp;

//	return result
	return prod;
}




////////////////////////////////////////////////////////////////////////
//	template instantiations for all current algebra types.
template class LocalSchurComplement<CPUAlgebra>;
template class LocalSchurComplement<CPUBlockAlgebra<3> >;
template class SchurComplementInverse<CPUAlgebra>;
template class SchurComplementInverse<CPUBlockAlgebra<3> >;
template class FETISolver<CPUAlgebra>;
template class FETISolver<CPUBlockAlgebra<3> >;
}// end of namespace
