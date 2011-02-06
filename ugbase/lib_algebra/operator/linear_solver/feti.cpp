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
//	used to inform root over primal connections.
template <class TValue>
struct PrimalConnection{
	PrimalConnection() {}
	PrimalConnection(int i1, int i2, TValue val) :
		ind1(i1), ind2(i2), value(val) {}
	int ind1;
	int ind2;
	TValue value;
};

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
	if(m_pOperator == NULL) // muesste man hier nicht pruefen, ob m_pDirichletMatrix gesetzt ist? -> @Ingo: Falsche Stelle:m_pDirichletMatrix ist die Matrix von m_DirichletOperator.
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
	uTmp = f;
	f.set(0.0);
	m_pFetiLayouts->vec_scale_append_on_dual(f, uTmp, 1.0);

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
	m_pCoarseProblemSolver(NULL),
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
	m_primalRootProc = 0;//pcl::GetOutputProcRank();

//	vector to store newly created root ids
//	please note that rootIDs may only be indexed with the agebra-index
//	of primal nodes. Content is not defined for other entries.
	std::vector<int> rootIDs;

//	Build layouts such that all processes can communicate their unknowns
//	to the primal Root Process
	int newVecSize = BuildOneToManyLayout(m_masterAllToOneLayout,
						 m_slaveAllToOneLayout, m_primalRootProc,
						 m_pFetiLayouts->get_primal_master_layout(),
						 m_pFetiLayouts->get_primal_slave_layout(),
						 m_allToOneProcessComm, &rootIDs);
/*
	UG_LOG("received root ids:");
	for(size_t i = 0; i < rootIDs.size(); ++i){
		UG_LOG(rootIDs[i] << " ");
	}
	UG_LOG(endl);
*/
//	We have to gather the rootIDs and the quantities of primal nodes
//	on each process of the feti-block in one array.
//	first we collect all local primal variables in one array
//	Collect all Primal indices on proc
	std::vector<IndexLayout::Element> vlocalPrimalIndex;
	CollectUniqueElements(vlocalPrimalIndex, m_slaveAllToOneLayout);


//	now build an array of local primal root ids
	std::vector<int> vlocalPrimalRootIDs(vlocalPrimalIndex.size());
	for(size_t i = 0; i < vlocalPrimalIndex.size(); ++i)
		vlocalPrimalRootIDs[i] = rootIDs[vlocalPrimalIndex[i]];

	pcl::ProcessCommunicator& localFetiBlockComm =
						m_pFetiLayouts->get_inner_process_communicator();

//	rootIDs of primal variables in the local feti block
	std::vector<int>	vPrimalQuantities;
//	vector that holds quantities of primal variables on each process
//	of the local feti block.
	std::vector<int> vPrimalRootIDs;

	localFetiBlockComm.allgatherv(vPrimalRootIDs, vlocalPrimalRootIDs,
								 &vPrimalQuantities);


//	log num primal quantities
	UG_LOG("primal root ids: ");
	for(size_t i = 0; i < vPrimalRootIDs.size(); ++i){
		UG_LOG(vPrimalRootIDs[i] << " ");
	}
	UG_LOG(std::endl);
	UG_LOG("primal quantities: ");
	for(size_t i = 0; i < vPrimalQuantities.size(); ++i){
		UG_LOG(vPrimalQuantities[i] << " ");
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
//	processes will collect their local primal connections here.
	typedef PrimalConnection<typename vector_type::value_type> PrimalConnection;
	vector<PrimalConnection> localPrimalConnections;

	vector_type e;  e.resize(m_pMatrix->num_rows());
	// set layouts and communicators for first help vectors (and duplicate this
	// by copying, instead of setting it over and over again for the other help vectors)
	m_pFetiLayouts->vec_use_std_communication(e); // TODO: correct communication??
	e.set(0.0);

	vector_type e2; e2.resize(m_pMatrix->num_rows()); e2 = e;
	vector_type e3; e3.resize(m_pMatrix->num_rows()); e3 = e;
	vector_type e4; e4.resize(m_pMatrix->num_rows()); e4 = e;
	vector_type e5; e5.resize(m_pMatrix->num_rows()); e5 = e;
	vector_type e6; e6.resize(m_pMatrix->num_rows()); e6 = e;

	size_t primalCounter = 0;
	for(size_t procInFetiBlock = 0; procInFetiBlock < localFetiBlockComm.size();
			procInFetiBlock++)
	{
		for(size_t pqi = 0; pqi < (size_t)vPrimalQuantities[procInFetiBlock]; ++pqi)
		{
		// 	1. Create unity vector
		//////////////////////////

			int connectedRootID = vPrimalRootIDs[primalCounter++];

		//	reset identity vector to zero for all Primal unknowns
			e.set(0.0);

		//	set value of unity vector to one if on process and quantity, else 0
			if(pcl::GetProcRank() == localFetiBlockComm.get_proc_id(procInFetiBlock))
			{
				const IndexLayout::Element localPrimalIndex = vlocalPrimalIndex[pqi];

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
						" local problem to compute Schur complement"
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


			// \todo: compute connectedRootID, i.e. root id of primal variable
			//			that has been set to one above

		//	loop process local primal unknowns
			for(size_t pqj = 0; pqj < vlocalPrimalIndex.size(); ++pqj)
			{
				const IndexLayout::Element localPrimalIndex = vlocalPrimalIndex[pqj];

				typename vector_type::value_type& entry = e6[localPrimalIndex];

				int primalRootID = rootIDs[localPrimalIndex];

				localPrimalConnections.push_back(PrimalConnection(connectedRootID,
															primalRootID, entry));
			}
		}
	}

//	all processes send their connections to root
//todo: This could be improved, so that only processes which contain
//		a primal node are involved.
	vector<PrimalConnection> vPrimalConnections;//	only filled on root
	pcl::ProcessCommunicator commWorld;
	commWorld.gatherv(vPrimalConnections, localPrimalConnections, m_primalRootProc);

//	log vPrimalConnections
	UG_LOG("primal connections:\n");
	for(size_t i = 0; i < vPrimalConnections.size(); ++i){
		PrimalConnection& pc = vPrimalConnections[i];
		UG_LOG("  ind1: " << pc.ind1 << "    ind2: " << pc.ind2 << "    value: " << pc.value << endl)
	}
	UG_LOG("endl");

// TODO: Hier muss die Schurkomplement-Matrix - m_pRootSchurComplementMatrix _ noch befuellt werden!?

// \todo: Remove, but currently here until invert is fully implemented
	return true;

//	init sequential solver for coarse problem
	if(m_pCoarseProblemSolver != NULL)
		if(!m_pCoarseProblemSolver->init(m_RootSchurComplementOp))
		{
			UG_LOG("ERROR in 'SchurComplementInverse::init': Cannot init "
					"coarse problem Solver for Operator S_{Pi Pi}.\n");
			return false;
		}


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
	vector_type fTmp; fTmp.create(f.size());
	vector_type hTmp; hTmp.create(u.size());

//	1. Set values of rhs to zero on I and Pi
	// (a) Reset all values - but not of original f, we need it later!
	fTmp.set(0.0);

	write_debug(f, "SCI_f_1BeforeNeumann");

	// (b) Copy values on \Delta
	m_pFetiLayouts->vec_scale_append_on_dual(fTmp, f, 1.0);

//	2. Compute \f$\tilde{f}_{\Pi}^{(p)}\f$ by computing \f$h_{\{I \Delta\}}^{(p)}\f$:
	write_debug(f,    "SCI_f_2aBeforeNeumann");
	write_debug(fTmp, "SCI_fTmp_2aBeforeNeumann");
	write_debug(hTmp, "SCI_hTmp_2aBeforeNeumann");

	// use inner interfaces for solving
	m_pFetiLayouts->vec_use_inner_communication(fTmp);
	fTmp.set_storage_type(PST_ADDITIVE);
	m_pFetiLayouts->vec_use_inner_communication(hTmp);
	hTmp.set_storage_type(PST_CONSISTENT);

	// (a) invoke Neumann solver to get \f$h_{\{I \Delta\}}^{(p)}\f$
	if(!m_pNeumannSolver->apply_return_defect(hTmp, fTmp))
	{
		UG_LOG_ALL_PROCS("ERROR in 'SchurComplementInverse::apply': "
						 "Could not solve Neumann problem (step 2.a) on Proc "
							<< pcl::GetProcRank() << ".\n");
		return false;
	}

	write_debug(hTmp, "SCI_hTmp_2a");

	// (b) apply matrix to \f$[h_{\{I \Delta\}}^{(p)}, 0]^T\f$ - multiply with full matrix
	if(!m_pMatrix->apply(fTmp, hTmp))
	{
		UG_LOG_ALL_PROCS("ERROR in 'SchurComplementInverse::apply': "
						 "Could not apply full matrix (step 2.b) on "
						 "Proc " << pcl::GetProcRank() << ".\n");
		return false;
	}

	write_debug(fTmp, "SCI_fTmp_2aAfterNeumann");

	// (c) Scale result by -1
	// (d) Set values to zero on I and Delta - excluding Pi
	hTmp = fTmp;
	fTmp.set(0.0);
	m_pFetiLayouts->vec_scale_append_on_primal(fTmp, hTmp, -1);

	write_debug(fTmp, "SCI_fTmp_2d");

//	3. Since \f$\tilde{f}_{\Pi}\f$ is saved additively, gather it to one process (root)
//     where it is then consistent.
	// TODO: Gather \f$\tilde{f}_{\Pi}\f$!

//	4. Solve \f$S_{\Pi \Pi} u_{\Pi} = \tilde{f}_{\Pi}\f$ on root
//     Toselli, p.~165, below eq.~(6.64): this is a ``local problem with Neumann bnd cnds at edges, zero Dirichlet bnd cnds''

	// TODO: Solve with 'SchurComplementInverse::m_pCoarseProblemSolver'!

//	5. Broadcast \f$u_{\Pi}\f$ to all Procs. \f$u_{\Pi}\f$ is consistently saved.
	// TODO: Broadcast solution!

	// Assumption: uTmp contains [0, 0, u_{\Pi}]^T
	vector_type uTmp; uTmp.create(f.size()); // den Hilfsvektor kann man sich vermutlich sparen ...
	// TODO: uTmp = ...
	// TODO: handle correctly
	uTmp.set(0.0);

//	6. Compute (via backward substitution)
	// \f$\hat{f}_{\{I \Delta\}}\f$ = *urspruengliches* f_{\{I \Delta\}} - Neumann-Matrix * [0, u_{\Pi}]^T
	m_pFetiLayouts->vec_use_inner_communication(hTmp);
	hTmp.set_storage_type(PST_ADDITIVE);
	m_pFetiLayouts->vec_use_inner_communication(uTmp);
	uTmp.set_storage_type(PST_CONSISTENT);
	if(!m_NeumannOperator.apply(hTmp, uTmp))
	{
		UG_LOG_ALL_PROCS("ERROR in 'SchurComplementInverse::apply': "
						 "Could not apply full matrix (step 6) on "
						 "Proc " << pcl::GetProcRank() << ".\n");
		return false;
	}

	write_debug(hTmp, "SCI_hTmp_6");

	VecScaleAdd(fTmp, 1.0, f, -1, hTmp);
	m_pFetiLayouts->vec_set_on_primal(fTmp, 0.0);

	write_debug(fTmp, "SCI_fTmp_6");
 
//	7. Solve for \f$u_{\{I \Delta\}}^{(p)}\f$
	m_pFetiLayouts->vec_use_inner_communication(fTmp);
	fTmp.set_storage_type(PST_ADDITIVE);
	m_pFetiLayouts->vec_use_inner_communication(uTmp);
	uTmp.set_storage_type(PST_CONSISTENT);
	if(!m_pNeumannSolver->apply_return_defect(uTmp, fTmp)) // solve with Neumann matrix!
	{
		UG_LOG_ALL_PROCS("ERROR in 'SchurComplementInverse::apply': "
						 "Could not solve Neumann problem (step 7) on Proc "
							<< pcl::GetProcRank() << ".\n");
		return false;
	}

	write_debug(uTmp, "SCI_uTmp_7");

//	ATTENTION: we skip this here, callers of this method must cancel values
//				if they need it. Otherwise we could not "backsolve"
//	8. Set values to zero on I and Pi - by excluding Delta
//	u.set(0.0);
//	m_pFetiLayouts->vec_scale_append_on_dual(u, uTmp, 1.0);
	u = uTmp;

	write_debug(u, "SCI_u_8");

//	we're done
	return true;
} /* end 'SchurComplementInverse::apply_return_defect()' */

template <typename TAlgebra>
bool SchurComplementInverse<TAlgebra>::
apply(vector_type& x, const vector_type& b)
{
	write_debug(b, "SCI_apply_b");

//	copy defect
	vector_type d; d.resize(b.size());
	d = b;

	write_debug(d, "SCI_apply_d");

//	solve on copy of defect
	return apply_return_defect(x, d);
} /* end 'SchurComplementInverse::apply()' */


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	FETISolver implementation
template <typename TAlgebra>
FETISolver<TAlgebra>::
FETISolver() :
	m_pOperator(NULL),
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
	m_pOperator = &A;

//	0. get matrix
	m_pMatrix = &m_pOperator->get_matrix();

//	check that DDInfo has been set
	if(m_pDDInfo == NULL)
	{
		UG_LOG("ERROR in FETISolver::init: DDInfo not set.\n");
		return false;
	}

	bool debugLayouts = (m_pDebugWriter==NULL) ? false : true;

//	1. create FETI Layouts
	m_fetiLayouts.create_layouts(m_pMatrix->get_master_layout(),
	                             m_pMatrix->get_slave_layout(),
	                             m_pMatrix->get_process_communicator(),
	                             m_pMatrix->num_rows(),
	                             *m_pDDInfo,
	                             debugLayouts);

//  ----- 2. CONFIGURE LOCAL SCHUR COMPLEMENT  ----- //

//	2.1 set layouts in LocalSchurComplement
	m_LocalSchurComplement.set_feti_layouts(m_fetiLayouts);

//	2.2 init Dirichlet system and solver
//  check that dirichlet solver has been set
	if(m_pDirichletSolver == NULL)
	{
		UG_LOG("ERROR in FETISolver::init: No dirichlet solver set "
				" for inversion of A_{II} in Local Schur complement.\n");
		return false;
	}

//	set dirichlet solver for local Schur complement
	m_LocalSchurComplement.set_dirichlet_solver(*m_pDirichletSolver);

//	set operator in local Schur complement
	m_LocalSchurComplement.set_matrix(*m_pOperator);

//	2.3 init local Schur complement
	if(m_LocalSchurComplement.init() != true)
	{
		UG_LOG("ERROR in FETISolver::init: Can not init local Schur "
				"complement.\n");
		bSuccess = false;
	}

//	2.4 check all procs
	if(!pcl::AllProcsTrue(bSuccess))
	{
		UG_LOG("ERROR in FETISolver::init: Some processes could not init"
				" local Schur complement.\n");
		return false;
	}

//  ----- 3. CONFIGURE SCHUR COMPLEMENT INVERSE  ----- //

//	3.1 set layouts in SchurComplementInverse
	m_SchurComplementInverse.set_feti_layouts(m_fetiLayouts);

//	3.2 init Neumann system and solver
//	check that neumann solver has been set
	if(m_pNeumannSolver == NULL)
	{
		UG_LOG("ERROR in FETISolver::init: No neumann solver set "
				" for inversion of A_{I,Delta}{I,Delta} in SchurComplementInverse.\n");
		return false;
	}

//	set neumann solver in SchurComplementInverse
	m_SchurComplementInverse.set_neumann_solver(*m_pNeumannSolver);

//  3.3 init coarse problem solver used in SchurComplementInverse
//	check that coarse problem solver has been set
	if(m_pCoarseProblemSolver == NULL)
	{
		UG_LOG("ERROR in FETISolver::init: No coarse problem solver set "
				" for solving S_{Pi Pi} u_{Pi} = tilde{f}_{Pi}.\n");
		return false;
	}

//	set coarse problem solver in SchurComplementInverse
	m_SchurComplementInverse.set_coarse_problem_solver(*m_pCoarseProblemSolver);

//	3.4 init SchurComplementInverse (operator - given as parameter here - is also set thereby)
	if(m_SchurComplementInverse.init(*m_pOperator) != true)
	{
		UG_LOG("ERROR in FETISolver::init: Can not init Schur "
				"complement inverse.\n");
		bSuccess = false;
	}

//	3.5 check all procs
	if(!pcl::AllProcsTrue(bSuccess))
	{
		UG_LOG("ERROR in FETISolver::init: Some process could not init"
				" Schur complement inverse.\n");
		return false;
	}

//  ----- 4. INIT DUAL DIRICHLET SOLVER  ----- //

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
	                                            m_fetiLayouts.get_primal_master_layout(),
	                                            m_fetiLayouts.get_primal_slave_layout());

//	init sequential solver for Dirichlet problem
	if(m_pDualDirichletSolver != NULL)
		if(!m_pDualDirichletSolver->init(m_DualDirichletOperator))
		{
			UG_LOG("ERROR in 'FETISolver::init': Cannot init "
					"Dual Dirichlet Solver for Operator A.\n");
			return false;
		}

//	Debug output of matrices
	if(m_pDebugWriter != NULL)
	{
		m_pDebugWriter->write_matrix(m_DualDirichletOperator.get_matrix(),
									 "FetiDualDirichletMatrix");

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
//	use storage for a whole vector and only use its Dual entries.

//	lagrange multiplier
	vector_type lambda; lambda.create(u.size());
	m_fetiLayouts.vec_use_std_communication(lambda);

//	residuum
	vector_type r; r.create(u.size());
	m_fetiLayouts.vec_use_std_communication(r);

//	search direction
	vector_type p; p.create(u.size());
	m_fetiLayouts.vec_use_std_communication(p);

//	preconditioned residuum
	vector_type z; z.create(u.size());

//	help vector to compute t = F*p
	vector_type t; t.create(u.size());
	m_fetiLayouts.vec_use_std_communication(t);

//	help vector to compute \tilde{f}_\Delta
	vector_type tildeF; tildeF.create(u.size());
	m_fetiLayouts.vec_use_std_communication(tildeF);


//	reset start value of lagrange multiplier
	lambda.set(0.0);

//	reset iteration count
	m_iterCnt = 0;

//	The idea of the Feti solver is to perform a cg-method on the Dual unknowns.
//	Therefore, the right-hand side has to be adapted, such that f_Dual fits to
//	the Schur complement system on the Dual unknowns. i.e.
//	\tilde{f}_{\Delta} = f_{\Delta}
//					   - A_{\Delta \{I \Pi\}} (A_{\{I \Pi\} \{I \Pi\}})^{-1} f_{\{I \Pi\}}
//
//	Please note: To compute this, we have to make the matrix consistent in the
//				 Primal unknowns. Then, the rhs can be computed on each process
//				 individually.

//	make f consistent on Pi
	AdditiveToConsistent(&f, m_fetiLayouts.get_primal_master_layout(),
	                     	 m_fetiLayouts.get_primal_slave_layout());

//	set inner layouts
	m_fetiLayouts.vec_use_inner_communication(tildeF);
	m_fetiLayouts.vec_use_inner_communication(f);

	write_debug(f, "ResiduumF");

//	compute \tilde{f}_{\Delta}
	if(!compute_tilde_f(tildeF, f))
	{
		UG_LOG("ERROR in 'FETISolver::apply': Unable "
			   "to build tilde f. Aborting.\n"); return false;
	}

	write_debug(tildeF, "ResiduumTildeF");

// 	Build start residuum:  r = r0 := d - F*lambda.
//	This is done in three steps.

//	(a) Compute d:= B S_{\Delta \Delta}^{-1} \tilde{f}_{\Delta}
//		and set r0 = d (prelimilarly, -F \lambda added later)
	tildeF.set_storage_type(PST_ADDITIVE);
	if(!compute_d(r, tildeF))
	{
		UG_LOG("ERROR in 'FETISolver::apply': "
			   "Cannot compute rhs 'd'. Aborting.\n");
		return false;
	}

	write_debug(r, "ResiduumD");

// 	(b) Build t = F*lambda (t is additive afterwards)
	UG_LOG(" ********** 'FETISOLVER::apply_return_defect()': Before 'apply_F()' ********** \n")
	if(!apply_F(t, lambda))
	{
		UG_LOG("ERROR in 'FETISolver::apply': Unable "
			   "to build t = F*p. Aborting.\n"); return false;
	}

// (c) Subtract values on \Delta, r0 = r0 - t
	m_fetiLayouts.vec_scale_append_on_dual(r, t, -1.0);

//	prepare appearance of conv check
	prepare_conv_check();

	write_debug(r, "Residuum");

//	compute and set start defect
	m_pConvCheck->start_defect(m_fetiLayouts.vec_norm_on_dual(r));

// 	Precondition the start defect: apply z = M^-1 * r
	UG_LOG(" ********** 'FETISOLVER::apply_return_defect()': Before 'apply_M_inverse()' ********** \n")
	if (!apply_M_inverse(z, r))
	{
		UG_LOG("ERROR in 'FETISolver::apply': "
			   "Cannot apply preconditioner. Aborting.\n");
		return false;
	}

//	start values
	number rho, rho_new, beta, alpha, alpha_denominator;
	rho = rho_new = beta = alpha = alpha_denominator = 0.0;

// 	start search direction
	p = z;

// 	start rho
	rho = m_fetiLayouts.vec_prod_on_dual(z, r);

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
		if (!apply_M_inverse(z, r))
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

	//	TODO:
	// 2. With \f$\lambda\f$ found, back solve for \f$u_{\Delta}\f$:
	//    \f$u_{\Delta} = {\tilde{S}_{\Delta \Delta}}^{-1} ({\tilde{f}_{\Delta}} - B_{\Delta}^T \lambda).\f$
	// 3. Make sure that all parts of the solution (\f$u_{I}, u_{\Delta}, $u_{\Pi}\f$)
	//    is assembled into the global solution vector


	//	"back solve" (cf. A. Toselli, O. Widlund: "Domain Decomposition Methods -
	//	Algorithms and Theory", cap. 6.4, p.164, l. 2) (03022011av)

//	reset t = 0.0
	t.set(0.0);

// compute t = B^T * lambda
	ComputeDifferenceOnDeltaTransposed(t, lambda,
	                                   m_fetiLayouts.get_dual_master_layout(),
	                                   m_fetiLayouts.get_dual_slave_layout(),
	                                   m_fetiLayouts.get_dual_nbr_slave_layout());

//	compute f = f - B^T * lambda
	m_fetiLayouts.vec_scale_append_on_dual(f, t, -1.0);

	write_debug(t, "FETI_t_Before_Sol");

//	Solve: A u = f
	if(!m_SchurComplementInverse.apply(u, f))
	{
		UG_LOG("ERROR in FETISolver::apply: Cannot back solve.\n");
		return false;
	}

	write_debug(u, "FETI_Sol");

	return m_pConvCheck->post();
} /* end 'FETISolver::apply_return_defect()' */

template <typename TAlgebra>
bool FETISolver<TAlgebra>::
apply_F(vector_type& f, const vector_type& v)
{
	//	Help vector
	vector_type fTmp; fTmp.create(v.size());
	//	0. Reset values of f, fTmp
	f.set(0.0); //fTmp.set(0.0);
	fTmp = f;
	// not nec. if tmp vector is copied:
	//fTmp.set_slave_layout(f.get_slave_layout());
	//fTmp.set_master_layout(f.get_master_layout());
	//fTmp.set_process_communicator(f.get_process_communicator());


	//	1. Apply transposed jump operator: f = B_{\Delta}^T * v_{\Delta}:
	ComputeDifferenceOnDeltaTransposed(f, v, m_fetiLayouts.get_dual_master_layout(),
	                                   	   	 m_fetiLayouts.get_dual_slave_layout(),
	                                   	   	 m_fetiLayouts.get_dual_nbr_slave_layout());

	//  2. Apply SchurComplementInverse to f
	m_SchurComplementInverse.apply(fTmp, f);

	//	3. Apply jump operator to get the final 'f'
	ComputeDifferenceOnDelta(f, fTmp, m_fetiLayouts.get_dual_master_layout(),
	                         	 	  m_fetiLayouts.get_dual_slave_layout(),
	                         	 	  m_fetiLayouts.get_dual_nbr_master_layout(),
	                         	 	  m_fetiLayouts.get_dual_nbr_slave_layout());

	//	we're done
	return true;
}


template <typename TAlgebra>
bool FETISolver<TAlgebra>::
compute_d(vector_type& d, const vector_type& f)
{
//	On entry, the vector f is filled with values for the dual unknowns. We make
//	no assumption on the values in the others (neglected unknowns).

//	create a help vector
	vector_type dTmp; dTmp.create(f.size());

//	0. Reset values of d, dTmp
	d.set(0.0); dTmp = d;

	write_debug(f, "ComputeD_f");

//  1. Apply SchurComplementInverse to 'f'
	if(!m_SchurComplementInverse.apply(dTmp, f))
	{
		UG_LOG("In 'FETISolver::compute_d': Could not apply Schur"
				" complement inverse.\n");
		return false;
	}

	write_debug(dTmp, "ComputeD_dTmp");

//	2. Apply jump operator to get the final 'd'
	ComputeDifferenceOnDelta(d, dTmp, m_fetiLayouts.get_dual_master_layout(),
	                         	 	  m_fetiLayouts.get_dual_slave_layout(),
	                         	 	  m_fetiLayouts.get_dual_nbr_master_layout(),
	                         	 	  m_fetiLayouts.get_dual_nbr_slave_layout());

	write_debug(d, "ComputeD_d");

//	we're done
	return true;
}

template <typename TAlgebra>
bool FETISolver<TAlgebra>::
compute_tilde_f(vector_type& tildeF, const vector_type& f)
{
//	Help vector
	vector_type fTmp; fTmp.create(f.size());
	vector_type hTmp; hTmp.create(f.size());
	fTmp = f;

//  set dirichlet zero values on f_Delta
	m_fetiLayouts.vec_set_on_dual(fTmp, 0.0);

	write_debug(fTmp, "ComputeTilde_fTmp");

//	solve system A_{\{I \PI\}\{I \PI\}}
	fTmp.set_storage_type(PST_ADDITIVE);
	hTmp.set_storage_type(PST_CONSISTENT);
	m_fetiLayouts.vec_use_inner_communication(hTmp);
	m_fetiLayouts.vec_use_inner_communication(fTmp);
	if(!m_pDualDirichletSolver->apply_return_defect(hTmp, fTmp))
	{
		UG_LOG("In 'FETISolver::compute_tilde_f': Could not apply inverse"
				" of Dual Dirichlet Matrix.\n");
		return false;
	}

	write_debug(hTmp, "ComputeTilde_hTmp");

//	multiply by A{\Delta {I \Pi}}
	if(!m_pOperator->apply(fTmp, hTmp))
	{
		UG_LOG("In 'FETISolver::compute_tilde_f': Could not apply Matrix.\n");
		return false;
	}

	write_debug(fTmp, "ComputeTilde_fTmpEnd");

//	build f := f - A*A^{-1}*f on dual
	m_fetiLayouts.vec_scale_add_on_dual(tildeF, 1.0, f, -1.0, fTmp);

	write_debug(tildeF, "ComputeTilde_tildeF");

//	we're done
	return true;
}

template <typename TAlgebra>
bool FETISolver<TAlgebra>::
apply_M_inverse(vector_type& z, const vector_type& r)
{
//	Help vector
	vector_type zTmp; zTmp.create(r.size()); zTmp = z;

//	0. Reset values of z, zTmp
	z.set(0.0); zTmp.set(0.0);

//	1. Apply scaling: z := D_{\Delta}^{(i)} * r
	apply_scaling_matrix(z, r); // maybe restrict to layout

//  2. Apply transposed jump operator: zTmp := B_{\Delta}^T * z
	ComputeDifferenceOnDeltaTransposed(zTmp, z, m_fetiLayouts.get_dual_master_layout(),
	                                   	   	    m_fetiLayouts.get_dual_slave_layout(),
	                                   	   	    m_fetiLayouts.get_dual_nbr_slave_layout());

//	3. Apply local Schur complement: z := S_{\Delta}^{(i)} * zTmp
//	3.1. let vectors use communication within feti subdomain
	m_fetiLayouts.vec_use_inner_communication(z);
	m_fetiLayouts.vec_use_inner_communication(zTmp);
//	3.2. set correct parallel storage type
	z.set_storage_type(PST_ADDITIVE);
	zTmp.set_storage_type(PST_CONSISTENT);
//	3.3. solve
	if(!m_LocalSchurComplement.apply(z, zTmp))
	{
		UG_LOG("ERROR in FETISolver::apply_M_inverse: Could not apply"
				" local Schur complement. \n");
		return false;
	}

//  4. Apply jump operator:  zTmp :=  B_{\Delta} * z
	ComputeDifferenceOnDelta(zTmp, z, m_fetiLayouts.get_dual_master_layout(),
	                         	 	  m_fetiLayouts.get_dual_slave_layout(),
	                         	 	  m_fetiLayouts.get_dual_nbr_master_layout(),
	                         	 	  m_fetiLayouts.get_dual_nbr_slave_layout());

//	5. Apply scaling: z := D_{\Delta}^{(i)} * zTmp to get the final 'z'
	apply_scaling_matrix(z, zTmp); // maybe restrict to layout

//	we're done
	return true;
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



/*	Temporaerer code, der in SchurComplementInverse eingebaut werden muss.
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
