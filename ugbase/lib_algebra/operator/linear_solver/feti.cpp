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
	m_pFetiLayouts(NULL),
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

//	check Feti layouts have been set
	if(m_pFetiLayouts == NULL)
	{
		UG_LOG("ERROR in 'LocalSchurComplement::init': FETI "
				" layouts not set.\n");
		return false;
	}

//	save matrix from which we build the Schur complement
	m_pMatrix = &m_pOperator->get_matrix();

//	get matrix from dirichlet operator
	m_pDirichletMatrix = &m_DirichletOperator.get_matrix();

//	Copy Matrix for Dirichlet Problem
	*m_pDirichletMatrix = *m_pMatrix;

//	Set Dirichlet values on Pi
	m_pFetiLayouts->mat_set_dirichlet_on_primal(*m_pDirichletMatrix);

//	Set Dirichlet values on Delta
	m_pFetiLayouts->mat_set_dirichlet_on_dual(*m_pDirichletMatrix);

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
	if(m_pOperator == NULL) // muesste man hier nicht pruefen, ob m_pDirichletMatrix gesetzt ist?
	{
		UG_LOG("ERROR: In 'LocalSchurComplement::apply': "
						"Matrix A not set.\n");
		return false;
	}

//	check Dirichlet solver
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
	m_pFetiLayouts->vec_scale_append_on_dual(uTmp, u, 1.0);

//	2. Compute rhs f_{I} = A_{I \Delta} u_{\Delta}
	if(!m_DirichletOperator.apply(f, uTmp))
	{
		UG_LOG_ALL_PROCS("ERROR in 'LocalSchurComplement::apply': "
						 "Could not compute Rhs for Dirichlet problem (step 2) on "
						 "Proc " << pcl::GetProcRank() << ".\n");
		return false;
	}
	// set values to zero on \Delta
	m_pFetiLayouts->vec_set_on_dual(f, 0.0);

//	3. Invert on inner unknowns u_{I} = A_{II}^{-1} f_{I}
	// (a) use the inner-FETI-block layouts
	//uTmp.use_layout(0);
	//f.use_layout(0);
	//m_pDirichletMatrix->use_layout(0);

	// (b) invoke Dirichlet solver
	if(!m_pDirichletSolver->apply_return_defect(uTmp, f))
	{
		UG_LOG_ALL_PROCS("ERROR in 'LocalSchurComplement::apply': "
						 "Could not solve Dirichlet problem (step 3.b) on Proc "
							<< pcl::GetProcRank() << ".\n");
		return false;
	}

//	4. Compute result vector
	// (a) Scale u_{I} by -1
	uTmp *= -1.0;

	// (b) Add u_{\Delta} on \Delta
	m_pFetiLayouts->vec_scale_append_on_dual(uTmp, u, 1.0);

	// (c) Multiply with full matrix
	if(!m_pOperator->apply(f, uTmp))
	{
		UG_LOG_ALL_PROCS("ERROR in 'LocalSchurComplement::apply': "
						 "Could not apply full matrix (step 4.c) on "
						 "Proc " << pcl::GetProcRank() << ".\n");
		return false;
	}

//	5. Reset all values for I, \Pi
	m_pFetiLayouts->vec_set_excl_dual(f, 0.0);

	// todo: this was written here in addition. IS IT NECESSARY/SENSEFUL -latuernich nicht! (Wer macht denn sowas??? ;-) ==> CAN BE DELETED!
//	VecSetOnLayout(&f, 0.0, m_pFetiLayouts->get_dual_slave_layout());
//	VecSetOnLayout(&f, 0.0, m_pFetiLayouts->get_dual_master_layout());

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
	m_pFetiLayouts(NULL),
	m_pNeumannMatrix(NULL),
	m_pNeumannSolver(NULL),
	m_primalRootProc(-1),
	m_pRootSchurComplementMatrix(NULL),
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
	if(m_pFetiLayouts == NULL)
	{
		UG_LOG_ALL_PROCS("ERROR in 'SchurComplementInverse::init':"
				" Feti Layouts not set "
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
	m_pFetiLayouts->mat_set_dirichlet_on_primal(*m_pNeumannMatrix);

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
	std::vector<int> rootIDs;//	-1 for all but for primal entries.

//	Build layouts such that all processes can communicated their unknowns
//	to the primal Root Process
	int newVecSize = BuildOneToManyLayout(m_masterAllToOneLayout,
						 m_slaveAllToOneLayout, m_primalRootProc,
						 m_pFetiLayouts->get_primal_master_layout(),
						 m_pFetiLayouts->get_primal_slave_layout(),
						 m_allToOneProcessComm, &rootIDs);

//	We have to gather the quantities of primal nodes on each process
//	of the feti-block in one array.
//todo: collect unique elements (remove below)
	int numLocalPrimals = m_slaveAllToOneLayout.num_interface_elements();
	pcl::ProcessCommunicator& localFetiBlockComm = m_pFetiLayouts->get_inner_process_communicator();
	UG_LOG("num local primals: " << numLocalPrimals << std::endl);
	m_primalQuantities.resize(localFetiBlockComm.size());
	localFetiBlockComm.allgather(&numLocalPrimals, 1, PCL_DT_INT,
							&m_primalQuantities.front(), 1, PCL_DT_INT);

//	todo: communicate newMasterIDs to all feti-block processes.
//		-> primalRootIDs (holds all root ids of primal variables of this
//		feti process.

//	log num primal quantities
	UG_LOG("proc ids: ");
	for(size_t i = 0; i < m_primalQuantities.size(); ++i){
		UG_LOG(localFetiBlockComm.get_proc_id(i) << " ");
	}
	UG_LOG(std::endl);
	UG_LOG("primal quantities: ");
	for(size_t i = 0; i < m_primalQuantities.size(); ++i){
		UG_LOG(m_primalQuantities[i] << " ");
	}
	UG_LOG(std::endl);

//	build matrix on primalRoot
	if(pcl::GetProcRank() == m_primalRootProc)
	{
	//	get matrix
		m_pRootSchurComplementMatrix = &m_RootSchurComplementOp.get_matrix();

	//	create matrix of correct size
		m_pRootSchurComplementMatrix->create(newVecSize, newVecSize);

		std::cout << "On PrimalRoot: Creating proc local Schur Complement"
					" of size " << newVecSize <<"x"<<newVecSize << std::endl;
	}

//	\todo: Compute matrix

	struct PrimalConnection{
		int ind1;
		int ind2;
		number value;
	};

	vector_type e; e.resize(m_pMatrix->num_rows());
	vector_type e2; e2.resize(m_pMatrix->num_rows());
	vector_type e3; e3.resize(m_pMatrix->num_rows());
	vector_type e4; e4.resize(m_pMatrix->num_rows());
	vector_type e5; e5.resize(m_pMatrix->num_rows());
	vector_type e6; e6.resize(m_pMatrix->num_rows());

	std::vector<IndexLayout::Element> vlocalPrimalIndex;

//	Collect all Primal indices on proc
	CollectUniqueElements(vlocalPrimalIndex, m_slaveAllToOneLayout);

	for(size_t procInFetiBlock = 0; procInFetiBlock < localFetiBlockComm.size();
			procInFetiBlock++)
	{
		for(size_t pqi = 0; pqi < (size_t)m_primalQuantities[procInFetiBlock]; ++pqi)
		{
		// 	1. Create unity vector
		//////////////////////////

		//	reset identity vector to zero for all Primal unknowns
			e.set(0.0);

		//	set value of unity vector to one if on process and quantity, else 0
			const IndexLayout::Element localPrimalIndex = vlocalPrimalIndex[pqi];
			if(pcl::GetProcRank() == localFetiBlockComm.get_proc_id(procInFetiBlock))
			{
				e[localPrimalIndex] = 1.0;
			}

		//	at this point the vector e has been set to an identity vector

		// 	2. Apply first matrix
		//////////////////////////
		//	build e2 = A_{\{I \Delta\} \Pi} e
			m_A->apply(e2, e);


		//	3. solve I,\Delta subsystem problem by:
		//////////////////////////
		//	Solve: A_{\{I \Delta\}  \{I \Delta\} } e3 = e2

		//	(a1) Set zero dirichlet bnd conds for e_2_{\Pi}
			m_pFetiLayouts->vec_set_on_primal(e2, 0.0);

		//	(a2) Start with zero iterate (not obligatory)
			e3.set(0.0);

		//	(b) Solve dirichlet problem
			if(!m_pNeumannSolver->apply(e3, e2))
			{
				UG_LOG("ERROR in SchurComplementInverse::init: Could not solve"
						" local neumann problem to compute Schur complement"
						" w.r.t. primal unknowns.\n");
				return false;
			}

		// 	4. Apply third matrix
		//////////////////////////

		//	(a) Set e3 zero on \Pi. This is enforced by neumann solver

		//	(b) apply matrix: e4 = A e3
			m_A->apply(e4, e3);

		//	(c) set entries to zero on I, \Delta (not needed, therefore skipped)

		// 	5. Compute first term
		///////////////////////////

		//	(a) multiply unity vector with matrix
			m_A->apply(e5, e);

		//	(b) reset values to zero on I, \Delta (not needed, therefore skipped)

		// 	6. Add values
		///////////////////////////

		//	e6 = e5 - e4
			m_pFetiLayouts->vec_scale_add_on_primal(e6, 1.0, e5, -1.0, e4);


		// 	at this point, we have the contribution of S_ij^{p} in all primal
		//	unknowns i. Thus, we have to read it and send it to the root process

//	hier senden wir besser noch nicht. Hier sollte man erstmal in eine struktur
//	struct PrimalConnections{
//		int masterID1, masterID2;
//		number val; - oder Matrix::entry_type oder so
//	};
//	und das in einen lokalen Vektor pushen. Gesendet wird am besten erst,
//	wenn die äußere Schleife beendet ist.

			int primalRootID = rootIDs[localPrimalIndex];

//			int connectedRootID
/*

			for(primalRootIDs...)
			{
			//	get entry
				typename vector_type::value_type& entry =
						e[GetEntryOfIndexLayout(m_slaveAllToOneLayout, pqj)];
				//\todo: continue implementation
			}
*/

		}
	}

//	ich denke hier sollte man senden... oder verstehe ich da was falsch?
//	...

//	we're done
	return true;
} /* end 'SchurComplementInverse::init()' */

template <typename TAlgebra>
bool SchurComplementInverse<TAlgebra>::
apply_return_defect(vector_type& u, vector_type& f)
{
//	check that matrix has been set
	if(m_pNeumannMatrix == NULL)
	{
		UG_LOG("ERROR: In 'SchurComplementInverse::apply': "
						"Matrix A not set.\n");
		return false;
	}

//	check Neumann solver
	if(m_pNeumannSolver == NULL)
	{
		UG_LOG("ERROR: In 'SchurComplementInverse::apply':"
						" No sequential Neumann Solver set.\n");
		return false;
	}

//	Check parallel storage type of matrix
	if(!m_pNeumannMatrix->has_storage_type(PST_ADDITIVE))
	{
		UG_LOG("ERROR: In 'SchurComplementInverse::apply': "
						"Inadequate storage format of matrix.\n");
		return false;
	}

//	Check parallel storage type of vectors
	if (!u.has_storage_type(PST_CONSISTENT))
	{
		UG_LOG("ERROR: In 'SchurComplementInverse::apply': "
						"Inadequate storage format of Vector 'u' (should be consistent).\n");
		return false;
	}
	if(!f.has_storage_type(PST_ADDITIVE))
	{
		UG_LOG("ERROR: In 'SchurComplementInverse::apply': "
						"Inadequate storage format of Vector 'f' (should be additive).\n");
		return false;
	}

//	Help vector
//	\todo: it would be sufficient to copy only the layouts without copying the values
	vector_type fTmp; fTmp.create(f.size()); fTmp = f;
	vector_type hTmp; hTmp.create(u.size()); hTmp = u;


//	1. Set values of rhs to zero on I and Pi
	// (a) Reset all values - but not of original f, we need it later!
	fTmp.set(0.0);

	// (b) Copy values on \Delta
	m_pFetiLayouts->vec_scale_append_on_dual(fTmp, f, 1.0);

//	2. Compute \f$\tilde{f}_{\Pi}^{(p)}\f$ by computing \f$h_{\{I \Delta\}}^{(p)}\f$:

	// (a) invoke Neumann solver to get \f$h_{\{I \Delta\}}^{(p)}\f$
	if(!m_pNeumannSolver->apply_return_defect(hTmp, fTmp))
	{
		UG_LOG_ALL_PROCS("ERROR in 'SchurComplementInverse::apply': "
						 "Could not solve Neumann problem (step 2.a) on Proc "
							<< pcl::GetProcRank() << ".\n");
		return false;
	}
	// (b) apply matrix to \f$[h_{\{I \Delta\}}^{(p)}, 0]^T\f$ - multiply with full matrix
	if(!m_pMatrix->apply(fTmp, hTmp))
	{
		UG_LOG_ALL_PROCS("ERROR in 'SchurComplementInverse::apply': "
						 "Could not apply full matrix (step 2.b) on "
						 "Proc " << pcl::GetProcRank() << ".\n");
		return false;
	}

	// (c) Scale result by -1
	fTmp *= -1.0;
	// (d) Set values to zero on I and Delta - excluding Pi
	m_pFetiLayouts->vec_set_excl_primal(fTmp, 0.0);

//	3. Since \f$\tilde{f}_{\Pi}\f$ is saved additively, gather it to one process (root)
//     where it is then consistent.
	// TODO: Gather \f$\tilde{f}_{\Pi}\f$!

//	4. Solve \f$S_{\Pi \Pi} u_{\Pi} = \tilde{f}_{\Pi}\f$ on root
	// TODO: Solve! Hinweis von Andreas: F¸r das Loesen der Matrix auf root beliebigen Loeser verwenden.
	// Dazu Objekt von ILinearOperatorInverse von aussen aufsetzen, im Skript z.B. LU (in Zukunft z.B. HLib ...)

//	5. Broadcast \f$u_{\Pi}\f$ to all Procs. \f$u_{\Pi}\f$ is consistently saved.
	// TODO: Broadcast solution!

	// Assumption: uTmp contains [0, 0, u_{\Pi}]^T
	vector_type uTmp; uTmp.create(f.size()); // den Hilfsvektor kann man sich vermutlich sparen ...
	// TODO: uTmp = ...

//	6. Compute (via backward substitution)
	// \f$\hat{f}_{\{I \Delta\}}\f$ = *urspruengliches* f_{\{I \Delta\}} - Neumann-Matrix * [0, u_{\Pi}]^T
	if(!m_pNeumannMatrix->apply(hTmp, uTmp))
	{
		UG_LOG_ALL_PROCS("ERROR in 'SchurComplementInverse::apply': "
						 "Could not apply full matrix (step 6) on "
						 "Proc " << pcl::GetProcRank() << ".\n");
		return false;
	}
	// war total falsch: m_pFetiLayouts->vec_scale_add_on_dual(fTmp, 1.0, f, -1.0, hTmp);
	VecScaleAdd(fTmp, 1.0, f, -1, hTmp);
	m_pFetiLayouts->vec_set_on_primal(fTmp, 0.0);

 
//	7. Solve for \f$u_{\{I \Delta\}}^{(p)}\f$
	if(!m_pNeumannSolver->apply_return_defect(u, fTmp)) // solve with Neumann matrix!
	{
		UG_LOG_ALL_PROCS("ERROR in 'SchurComplementInverse::apply': "
						 "Could not solve Neumann problem (step 7) on Proc "
							<< pcl::GetProcRank() << ".\n");
		return false;
	}

//	8. Set values to zero on I and Pi - by excluding Delta
	m_pFetiLayouts->vec_set_excl_dual(u, 0.0);

//	we're done
	return true;
} /* end 'SchurComplementInverse::apply_return_defect()' */

template <typename TAlgebra>
bool SchurComplementInverse<TAlgebra>::
apply(vector_type& x, const vector_type& b)
{
	
//	copy defect
	vector_type d; d.resize(b.size());
	d = b;

//	solve on copy of defect
	return apply_return_defect(x, d);
} /* end 'SchurComplementInverse::apply()' */


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

//	check that DDInfo has been set
	if(m_pDDInfo == NULL)
	{
		UG_LOG("ERROR in FETISolver::init: DDInfo not set.\n");
		return false;
	}

//	create FETI Layouts
	m_fetiLayouts.create_layouts(m_pMatrix->get_master_layout(),
	                             m_pMatrix->get_slave_layout(),
	                             m_pMatrix->num_rows(),
	                             *m_pDDInfo);

//  ----- INIT DIRICHLET SOLVER  ----- //

//	check that dirichlet solver has been set
	if(m_pDirichletSolver == NULL)
	{
		UG_LOG("ERROR in FETISolver::init: No dirichlet solver set "
				" for inversion of A_{II} in Local Schur complement.\n");
		return false;
	}

//	set layouts
	m_LocalSchurComplement.set_feti_layouts(m_fetiLayouts);

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

//	set layouts
	m_SchurComplementInverse.set_feti_layouts(m_fetiLayouts);

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
	m_fetiLayouts.mat_set_dirichlet_on_dual(*m_pDualDirichletMatrix);

//	make diagonal of matrix consistent on Pi
	MatAdditiveToConsistentOnDiag<algebra_type>(m_pDualDirichletMatrix,
	                                            m_fetiLayouts.get_dual_master_layout(),
	                                            m_fetiLayouts.get_dual_slave_layout());

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
	lambda.set_slave_layout(u.get_slave_layout());
	lambda.set_master_layout(u.get_master_layout());
	lambda.set_process_communicator(u.get_process_communicator());

//	residuum
	vector_type r; r.create(u.size());
	r.set_slave_layout(u.get_slave_layout());
	r.set_master_layout(u.get_master_layout());
	r.set_process_communicator(u.get_process_communicator());

//	search direction
	vector_type p; p.create(u.size());
	p.set_slave_layout(u.get_slave_layout());
	p.set_master_layout(u.get_master_layout());
	p.set_process_communicator(u.get_process_communicator());

//	preconditioned residuum
	vector_type z; z.create(u.size());

//	help vector to compute t = F*p
	vector_type t; t.create(u.size());
	t.set_slave_layout(u.get_slave_layout());
	t.set_master_layout(u.get_master_layout());
	t.set_process_communicator(u.get_process_communicator());


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
	AdditiveToConsistent(&f, m_fetiLayouts.get_primal_master_layout(),
	                     	 m_fetiLayouts.get_primal_slave_layout());

//	set inner layouts
	m_fetiLayouts.vec_use_inner_communication(t);
	m_fetiLayouts.vec_use_inner_communication(f);

//	compute \tilde{f}_{\Delta}
	if(!compute_tilde_f(t, f)) return false;

// 	Build start residuum:  r = r0 := d - F*lambda.
//	This is done in three steps.

//	(a) Compute d:= B S_{\Delta \Delta}^{-1} \tilde{f}_{\Delta}
//		and set r0 = d -- only if lambda^0 = 0!
	if(!compute_d(r, f))
	{
		UG_LOG("ERROR in 'FETISolver::apply': "
			   "Cannot compute rhs 'd'. Aborting.\n");
		return false;
	}

// 	(b) Build t = F*lambda (t is additive afterwards)
	UG_LOG(" ********** 'FETISOLVER::apply_return_defect()': Before 'apply_F()' ********** \n")
	if(!apply_F(t, lambda))
	{
		UG_LOG("ERROR in 'FETISolver::apply': Unable "
			   "to build t = F*p. Aborting.\n"); return false;
	}

// (c) Subtract values on \Delta, r0 = r0 - t
	m_fetiLayouts.vec_scale_append_on_dual(r, t, -1.0);

// 	Precondition the start defect: apply z = M^-1 * r
	UG_LOG(" ********** 'FETISOLVER::apply_return_defect()': Before 'apply_M_inverse_with_identity_scaling()' ********** \n")
	if (!apply_M_inverse_with_identity_scaling(z, r))
	{
		UG_LOG("ERROR in 'FETISolver::apply': "
			   "Cannot apply preconditioner. Aborting.\n");
		return false;
	}

//	prepare appearance of conv check
	prepare_conv_check();

//	compute and set start defect
	m_pConvCheck->start_defect(m_fetiLayouts.vec_norm_on_dual(r));

//	start values
	number rho, rho_new, beta, alpha, alpha_denominator;
	rho = rho_new = beta = alpha = alpha_denominator = 0.0;

// 	start search direction
	p = z;

// 	start rho
	rho = m_fetiLayouts.vec_prod_on_dual(z, r);

//	reset iteration count
	m_iterCnt = 0;

// 	Iteration loop
	while(!m_pConvCheck->iteration_ended())
	{
	//	increase iteration count
		m_iterCnt++;
		UG_LOG(" ********** 'FETISOLVER::apply_return_defect()': iter " << m_iterCnt << " ********** \n")

	// 	Build t = F*p (t is additive afterwards)
		if(!apply_F(t, p))
		{
			UG_LOG("ERROR in 'FETISolver::apply': Unable "
						"to build t = F*p. Aborting.\n"); return false;
		}

	// 	Compute alpha
		alpha_denominator = m_fetiLayouts.vec_prod_on_dual(t, p);

		alpha = rho/alpha_denominator;

	// 	Update lambda := lambda + alpha*p
		m_fetiLayouts.vec_scale_add_on_dual(lambda, 1.0, lambda, alpha, p);

	// 	Update r := r - alpha*t
		m_fetiLayouts.vec_scale_add_on_dual(r, 1.0, r, -alpha, t);

	// 	Compute new defect
		m_pConvCheck->update_defect(m_fetiLayouts.vec_norm_on_dual(r));
		UG_LOG(" ********** 'FETISOLVER::apply_return_defect()': iter " << m_iterCnt << ": After m_pConvCheck->update(r) ********** \n")

	// 	Preconditioning: apply z = M^-1 * r
		if (!apply_M_inverse_with_identity_scaling(z, r))
		{
			UG_LOG("ERROR in 'FETISolver::apply': "
				   "Cannot apply preconditioner. Aborting.\n"); return false;
		}

	// 	new rho
		rho_new = m_fetiLayouts.vec_prod_on_dual(z, r);

	// 	new beta
		beta = rho_new/rho;

	// 	new direction p:= beta*p + z
		m_fetiLayouts.vec_scale_add_on_dual(p, beta, p, 1.0, z);

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



/*	Temporärer code, der in SchurComplementInverse eingebaut werden muss.
struct PrimalConnection{
	int ind1;
	int ind2;
	number value;
};

void tmp()
{
std::vector<PrimalConnection> localConns;

//...

//	communicate sizes first
ProcessCommunicator comm;
std::vector<int> allConnSizes(comm.size(), 0);

int localSize = (int) localConns.size();
comm.gather(&localSize, 1, PCL_DT_INT, &allConnSizes.front(),
			1, MPI_DT_INT, rootProc);

int totalNumConns = 0;
std::vector<int> offsets(comm.size(), 0);
//	since we're sending the data as byte, we have to adjust
//	sizes and offsets.
for(size_t i = 0; i < allConnSizes.size(); ++i){
	offsets[i] = totalNumConns * sizeof(PrimalConnection);
	totalNumConns += allConnSizes[i];
	allConnSizes[i] *= sizeof(PrimalConnection);
}

//	root now knows all sizes. We can now gather the connections on root.
std::vector<PrimalConnection> allConns(totalNumConns);
comm.gatherv(&localConns.front(),
			 sizeof(PrimalConnection) * localConns.size(),
			 PCL_DT_BYTE,
			 &allConns.front(),
			 &allConnSizes.front(),
			 &offsets.front(),
			 MPI_DT_INT,
			 rootProc);

}
*/
