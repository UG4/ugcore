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

namespace ug{

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	LocalSchurComplement implementation
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
		m_pDebugWriter->write_matrix(m_DirichletOperator.get_matrix(),
									 "FetiDirichletMatrix");
		m_pDebugWriter->write_matrix(m_pOperator->get_matrix(),
									 "FetiOriginalMatrix");
	}

UG_LOG_ALL_PROCS("'LocalSchurComplement::init': "
		  "********** After MatSetDirichletOnLayout on Proc " << pcl::GetProcRank() << ".(TMP)\n");

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
	VecScaleAddOnLayout(&uTmp, &u, 1.0, *m_pSlaveDualLayout);
	VecScaleAddOnLayout(&uTmp, &u, 1.0, *m_pMasterDualLayout);

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
	VecScaleAddOnLayout(&uTmp, &u, 1.0, *m_pSlaveDualLayout);
	VecScaleAddOnLayout(&uTmp, &u, 1.0, *m_pMasterDualLayout);

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
bool SchurComplementInverse<TAlgebra>::
init(ILinearOperator<vector_type, vector_type>& L)
{
//	success flag
	bool bSuccess = true;

	pcl::SynchronizeProcesses();
	UG_LOG_ALL_PROCS("Initializing SchurComplementInverse.\n");

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
	m_primalQuantities.resize(m_localFetiBlockComm.size());
	m_localFetiBlockComm.allgather(&numLocalPrimals, 1, PCL_DT_INT,
							&m_primalQuantities.front(), 1, PCL_DT_INT);

//	log num primal quantities
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
					" of size " << newVecSize <<"^2." << std::endl;
	}

//	we're done
	return true;
} /* end 'SchurComplementInverse::init()' */



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	FETISolver implementation
template <typename TAlgebra>
bool FETISolver<TAlgebra>::
init(IMatrixOperator<vector_type, vector_type, matrix_type>& A)
{
//	bool flag
	bool bSuccess;

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

	pcl::SynchronizeProcesses();
	UG_LOG(" ********** INIT LOCAL SCHUR COMPLEMENT ********** ... \n")

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

	pcl::SynchronizeProcesses();
	UG_LOG(" ********** INIT LOCAL SCHUR COMPLEMENT DONE ********** \n")

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

	pcl::SynchronizeProcesses();
	UG_LOG(" ********** INIT SCHUR COMPLEMENT INVERSE ********** ... \n")

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

	pcl::SynchronizeProcesses();
	UG_LOG(" ********** INIT SCHUR COMPLEMENT INVERSE DONE ********** \n")

//	we're done
	return true;
} /* end 'FETISolver::init()' */

template <typename TAlgebra>
bool FETISolver<TAlgebra>::
apply_return_defect(vector_type& lambda, vector_type& d)
{
	if(m_pConvCheck == NULL)
	{
		UG_LOG("ERROR: In 'FETISolver::apply': Convergence check not set.\n");
		return false;
	}

	#ifdef UG_PARALLEL
	if(!d.has_storage_type(PST_ADDITIVE) || !lambda.has_storage_type(PST_CONSISTENT))
		{
			UG_LOG("ERROR: In 'FETISolver::apply':Inadequate storage format of Vectors 'd' and 'lambda'.\n");
			return false;
		}
	#endif

// 	rename r as d (for convenience)
	vector_type& r = d;

// 	create help vectors (h will be consistent r) -- 'h'???? Wieso unten sizes unterschiedlicher vectors im 'create()'? Die sollten doch alle gleich lang sein ...
//	todo: 	It would be sufficient to copy only the pattern and
//			without initializing, but in parallel we have to copy communicators
	vector_type t; t.create(r.size());      t = r;
	vector_type z; z.create(lambda.size()); z = lambda;
	vector_type p; p.create(lambda.size()); p = lambda;

// 	Build residuum:  r := d - F*lambda
	//if(!m_A->apply_sub(r, lambda))
// 	(a) Build t = F*lambda (t is additive afterwards)
	if(!apply_F(t, lambda))
	{
		UG_LOG("ERROR in 'FETISolver::apply': Unable "
			   "to build t = F*p. Aborting.\n"); return false;
	}
// (b) Copy values on \Delta
	VecScaleAddOnLayout(&r, &d, -1.0, m_masterDualLayout);
	VecScaleAddOnLayout(&r, &d, -1.0, m_slaveDualLayout);

// 	Preconditioning: apply z = M^-1 * r
	if (!apply_M_inverse_with_identity_scaling(z, r))
	{
		UG_LOG("ERROR in 'FETISolver::apply': "
			   "Cannot apply preconditioner. Aborting.\n"); return false;
	}
	else z = r;

	// make z consistent
	if(!z.change_storage_type(PST_CONSISTENT))
	{
		UG_LOG("ERROR in 'FETISolver::apply': "
			   "Cannot convert z to consistent vector.\n"); return false;
	}

	prepare_conv_check();
	m_pConvCheck->start(r);

	number rho, rho_new, beta, alpha, alpha_denominator;
	rho = rho_new = beta = alpha = alpha_denominator = 0.0;

	// start search direction
	p = z;

	// start rho
	//rho = VecProd(z, r);
	VecProdOnLayout(rho, &z, &r, m_masterDualLayout);

	m_iterCnt = 0;
// 	Iteration loop
	while(!m_pConvCheck->iteration_ended())
	{
		m_iterCnt++;
	// 	Build t = F*p (t is additive afterwards)
		if(!apply_F(t, p))
		{
			UG_LOG("ERROR in 'FETISolver::apply': Unable "
						"to build t = F*p. Aborting.\n"); return false;
		}

	// 	Compute alpha
		//alpha_denominator = VecProd(t, p);
		VecProdOnLayout(alpha_denominator, &t, &p, m_masterDualLayout);

		alpha = rho/alpha_denominator;

	// 	Update lambda := lambda + alpha*p
		//VecScaleAdd(lambda, 1.0, lambda, alpha, p);
		VecScaleAddOnLayout(&lambda, 1.0, &lambda, alpha, &p, m_slaveDualLayout);
		VecScaleAddOnLayout(&lambda, 1.0, &lambda, alpha, &p, m_masterDualLayout);

	// 	Update r := r - alpha*t
		//VecScaleAdd(r, 1.0, r, -alpha, t);
		VecScaleAddOnLayout(&r, 1.0, &r, -alpha, &t, m_slaveDualLayout);
		VecScaleAddOnLayout(&r, 1.0, &r, -alpha, &t, m_masterDualLayout);

	// 	Check convergence
		m_pConvCheck->update(r);

	// 	Preconditioning: apply z = M^-1 * r
		if (!apply_M_inverse_with_identity_scaling(z, r))
		{
			UG_LOG("ERROR in 'FETISolver::apply': "
				   "Cannot apply preconditioner. Aborting.\n"); return false;
		}
		else z = r;

	// 	make z consistent
		if(!z.change_storage_type(PST_CONSISTENT))
		{
			UG_LOG("ERROR in 'FETISolver::apply': "
				   "Cannot convert z to consistent vector.\n"); return false;
		}

	// 	new rho
		//rho_new = VecProd(z, r);
		VecProdOnLayout(rho_new, &z, &r, m_masterDualLayout);

	// 	new beta
		beta = rho_new/rho;

	// 	new direction p:= beta*p + z
		//VecScaleAdd(p, beta, p, 1.0, z);
		VecScaleAddOnLayout(&p, beta, &p, 1.0, &z, m_slaveDualLayout);
		VecScaleAddOnLayout(&p, beta, &p, 1.0, &z, m_masterDualLayout);

	// 	update rho
		rho = rho_new;
	} /* end iteration loop */

	return m_pConvCheck->post();
} /* end 'FETISolver::apply_return_defect()' */





////////////////////////////////////////////////////////////////////////
//	template instantiations for all current algebra types.
template class LocalSchurComplement<CPUAlgebra>;
template class LocalSchurComplement<CPUBlockAlgebra<3> >;
template class SchurComplementInverse<CPUAlgebra>;
template class SchurComplementInverse<CPUBlockAlgebra<3> >;
template class FETISolver<CPUAlgebra>;
template class FETISolver<CPUBlockAlgebra<3> >;
}// end of namespace
