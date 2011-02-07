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
	m_pFetiLayouts->vec_use_inner_communication(f);
	f.set_storage_type(PST_ADDITIVE);
	m_pFetiLayouts->vec_use_inner_communication(uTmp);
	uTmp.set_storage_type(PST_CONSISTENT);

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
//	PrimalSubassembledMatrixInverse implementation
template <typename TAlgebra>
PrimalSubassembledMatrixInverse<TAlgebra>::
PrimalSubassembledMatrixInverse() :
	m_pOperator(NULL),
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
bool PrimalSubassembledMatrixInverse<TAlgebra>::
init(ILinearOperator<vector_type, vector_type>& L)
{
//	success flag
	bool bSuccess = true;

//	remember operator
	m_pOperator = dynamic_cast<IMatrixOperator<vector_type, vector_type, matrix_type>*>(&L);

//	check, that operator is correct
	if(m_pOperator == NULL)
	{
		UG_LOG_ALL_PROCS("ERROR in 'PrimalSubassembledMatrixInverse::init':"
				" Wrong type of operator passed for init.\n");
		bSuccess = false;
	}

//	check that Pi layouts have been set
	if(m_pFetiLayouts == NULL)
	{
		UG_LOG_ALL_PROCS("ERROR in 'PrimalSubassembledMatrixInverse::init':"
				" Feti Layouts not set "
				"on Proc " << pcl::GetProcRank() << ".\n");
		bSuccess = false;
	}

//	Check all procs
	if(!pcl::AllProcsTrue(bSuccess))
	{
		UG_LOG_ALL_PROCS("ERROR in 'PrimalSubassembledMatrixInverse::init':"
				" Some proc could not init Schur Complement inverse.\n");
		return false;
	}

//	save matrix from which we build the Schur complement
	m_pMatrix = &m_pOperator->get_matrix();

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
			UG_LOG("ERROR in 'PrimalSubassembledMatrixInverse::init': Cannot init "
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
//	please note that rootIDs may only be indexed with the algebra-index
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

//	processes will collect their local primal connections here.
	typedef PrimalConnection<typename vector_type::value_type> PrimalConnection;
	vector<PrimalConnection> localPrimalConnections;

//	create help vectors
	vector_type e; e.resize(m_pMatrix->num_rows());
	vector_type h1; h1.resize(m_pMatrix->num_rows());
	vector_type h2; h2.resize(m_pMatrix->num_rows());

//	set communication to inner subdomain comm.
	m_pFetiLayouts->vec_use_inner_communication(e);
	m_pFetiLayouts->vec_use_inner_communication(h1);
	m_pFetiLayouts->vec_use_inner_communication(h2);

//	Now within each feti subdomain, the primal unknowns are loop one after the
//	other. This is done like the following: Each process loops the number of
//	procs of the feti subdomain, and then for each subdomain the number of
//	primal unknowns (as stored in vPrimalQuantities) on this subdomain. That way
//	for every primal unknown on the feti subdomain, we can set its value to 1
//	while all the other values are zero. This gives us the unity vector. For
//	this unity vector an application of S_{Pi Pi} is computed such that we can
//	read then the values of (S_{Pi Pi})_{ij} in the i-th component of the result
//	vector. All those couplings are stored in the vector of connections called
//	localPrimalConnections and sent to the primalRootProc at the end of the
//	loop.
	size_t primalCounter = 0;
	for(size_t procInFetiBlock = 0; procInFetiBlock < localFetiBlockComm.size();
			procInFetiBlock++)
	{
		for(size_t pqi = 0; pqi < (size_t)vPrimalQuantities[procInFetiBlock]; ++pqi)
		{
		////////////////////////////
		// 	1. Create unity vector
		////////////////////////////

		//	remember the matrix id on root, of the primal unknown, that is set to one
			int unityRootID = vPrimalRootIDs[primalCounter++];

		//	reset identity vector to zero for all Primal unknowns
			e.set(0.0);

		//	set value of unity vector to one if on process and quantity, else 0
			if(pcl::GetProcRank() == localFetiBlockComm.get_proc_id(procInFetiBlock))
			{
				const IndexLayout::Element localPrimalIndex = vlocalPrimalIndex[pqi];

				e[localPrimalIndex] = 1.0;
			}

		//////////////////////////
		// 	2. Apply first matrix
		//////////////////////////
		//	build h1 = A_{\{I \Delta\} \Pi} e
			m_pOperator->apply(h1, e);

		////////////////////////////////////////////
		//	3. solve I,\Delta subsystem problem by:
		////////////////////////////////////////////
		//	Solve: A_{\{I \Delta\}  \{I \Delta\} } h2 = h1

		//	(a1) Set zero dirichlet bnd conds for e_2_{\Pi}
			m_pFetiLayouts->vec_set_on_primal(h1, 0.0);

		//	(a2) Start with zero iterate (not obligatory)
			h2.set(0.0);

		//	(b) Solve dirichlet problem
			if(!m_pNeumannSolver->apply(h2, h1))
			{
				UG_LOG("ERROR in PrimalSubassembledMatrixInverse::init: Could not solve"
						" local problem to compute Schur complement"
						" w.r.t. primal unknowns.\n");
				return false;
			}

		//////////////////////////
		// 	4. Apply third matrix
		//////////////////////////

		//	(a) Set h2 zero on \Pi. This is enforced by neumann solver

		//	(b) apply matrix: h1 = A h2
			m_pOperator->apply(h1, h2);

		//	(c) set entries to zero on I, \Delta (not needed, therefore skipped)

		///////////////////////////
		// 	5. Compute first term
		///////////////////////////

		//	(a) multiply unity vector with matrix
			m_pOperator->apply(h2, e);

		//	(b) reset values to zero on I, \Delta (not needed, therefore skipped)

		///////////////////////////
		// 	6. Add values
		///////////////////////////

		//	e = h2 - h1
			m_pFetiLayouts->vec_scale_add_on_primal(e, 1.0, h2, -1.0, h1);


		// 	at this point, we have the contribution of S_ij^{p} in all primal
		//	unknowns i. Thus, we have to read it and send it to the root process

		//	loop process local primal unknowns
			for(size_t pqj = 0; pqj < vlocalPrimalIndex.size(); ++pqj)
			{
				const IndexLayout::Element localPrimalIndex = vlocalPrimalIndex[pqj];

			//	read coupling
				typename vector_type::value_type& entry = e[localPrimalIndex];

			//	read root index
				int primalRootID = rootIDs[localPrimalIndex];

			//  remember coupling
				localPrimalConnections.push_back(PrimalConnection(primalRootID,
				                                                  unityRootID, entry));
			}
		}
	}

//	all processes send their connections to root
//todo: This could be improved, so that only processes which contain
//		a primal node are involved.
	vector<PrimalConnection> vPrimalConnections;//	only filled on root
	pcl::ProcessCommunicator commWorld;
	commWorld.gatherv(vPrimalConnections, localPrimalConnections, m_primalRootProc);

//	build matrix on primalRoot
	if(pcl::GetProcRank() == m_primalRootProc)
	{
	//	get matrix
		m_pRootSchurComplementMatrix = &m_RootSchurComplementOp.get_matrix();

	//	check matrix
		if(m_pRootSchurComplementMatrix == NULL)
		{
			UG_LOG("ERROR in 'PrimalSubassembledMatrixInverse::init': No matrix in"
					"Root Schur Complement Operator.\n");
			return false;
		}

	//	reference for convinience
		matrix_type& mat = *m_pRootSchurComplementMatrix;

	//	create matrix of correct size
		mat.resize(newVecSize, newVecSize);

	//	info output
		std::cout << "On PrimalRoot: Creating proc local Schur Complement"
					" of size " << newVecSize <<"x"<<newVecSize << std::endl;

	//	copy received values into matrix
		mat.set(0.0);
		std::cout << "Writing primal connections:" << std::endl;
		for(size_t i = 0; i < vPrimalConnections.size(); ++i)
		{
		//	get sent connection
			PrimalConnection& pc = vPrimalConnections[i];
			std::cout << "  ind1: " << pc.ind1 << "    ind2: " << pc.ind2 << "    value: " << pc.value << std::endl;

		//	get corresponding block
			typename matrix_type::value_type& block = mat(pc.ind1, pc.ind2);

		//	loop block components
			for(size_t beta = 0; beta < (size_t) GetCols(block); ++beta)
			{
				BlockRef(block, beta, beta) += BlockRef(pc.value, beta);
			}
		}

	//	init sequential solver for coarse problem
		if(newVecSize > 0)
			if(m_pCoarseProblemSolver != NULL)
			{
				if(!m_pCoarseProblemSolver->init(m_RootSchurComplementOp))
				{
					UG_LOG("ERROR in 'PrimalSubassembledMatrixInverse::init': Cannot init "
							"coarse problem Solver for Operator S_{Pi Pi}.\n");
					return false;
				}
			}
			else
			{
				UG_LOG("ERROR in 'PrimalSubassembledMatrixInverse::init': S_{Pi Pi is} "
						" needs to be inverted, but no CoarseSolver given.\n");
				return false;
			}

	//	Debug output of matrix (only used to flag if debug is wanted)
		if(m_pDebugWriter != NULL)
		{
		//	use own debug writer
			AlgebraDebugWriter<algebra_type, 2> debugWriter;

		//	vector of positions
			std::vector<MathVector<2> > pos(newVecSize);
			if(newVecSize == 9)
			{
				pos[0] = MathVector<2>(-0.5, -0.5);
				pos[1] = MathVector<2>(-0.0, -0.5);
				pos[2] = MathVector<2>(-0.5, -0.0);
				pos[3] = MathVector<2>(0.0, 0.0);
				pos[4] = MathVector<2>(0.5, -0.5);
				pos[5] = MathVector<2>(0.5, 0.0);
				pos[6] = MathVector<2>(-0.5, 0.5);
				pos[7] = MathVector<2>(0.0, 0.5);
				pos[8] = MathVector<2>(0.5, 0.5);
			}
		//	set positions
			debugWriter.set_positions(&pos[0], newVecSize);

			debugWriter.write_matrix(m_RootSchurComplementOp.get_matrix(),
										 "RootSchurComplementMatrix");
		}
	}

//	we're done
	return true;
}
/* end 'PrimalSubassembledMatrixInverse::init()' */

template <typename TAlgebra>
bool PrimalSubassembledMatrixInverse<TAlgebra>::
apply_return_defect(vector_type& u, vector_type& f)
{
//	check that matrix has been set
	if(m_pNeumannMatrix == NULL)
	{
		UG_LOG("ERROR: In 'PrimalSubassembledMatrixInverse::apply': "
						"Matrix A not set.\n");
		return false;
	}

//	check Neumann solver
	if(m_pNeumannSolver == NULL)
	{
		UG_LOG("ERROR: In 'PrimalSubassembledMatrixInverse::apply':"
						" No sequential Neumann Solver set.\n");
		return false;
	}

//	Check parallel storage type of matrix
	if(!m_pNeumannMatrix->has_storage_type(PST_ADDITIVE))
	{
		UG_LOG("ERROR: In 'PrimalSubassembledMatrixInverse::apply': "
						"Inadequate storage format of matrix.\n");
		return false;
	}

//	Check parallel storage type of vectors
	if (!u.has_storage_type(PST_CONSISTENT))
	{
		UG_LOG("ERROR: In 'PrimalSubassembledMatrixInverse::apply': "
						"Inadequate storage format of Vector 'u' (should be consistent).\n");
		return false;
	}
	if(!f.has_storage_type(PST_ADDITIVE))
	{
		UG_LOG("ERROR: In 'PrimalSubassembledMatrixInverse::apply': "
						"Inadequate storage format of Vector 'f' (should be additive).\n");
		return false;
	}

//	success flag
	bool bSuccess = true;

//	0. create help vector
	vector_type h; h.create(u.size());

//	1. Set values of rhs to zero on Pi
	// (a) Copy values
	h = f;

	// (b) Reset values on Pi
	m_pFetiLayouts->vec_set_on_primal(h, 0.0);

//	2. Compute \f$\tilde{f}_{\Pi}^{(p)}\f$ by computing \f$h_{\{I \Delta\}}^{(p)}\f$:
	write_debug(h, "PSMI_h_2_BeforeNeumann");

	// use inner interfaces for solving
	m_pFetiLayouts->vec_use_inner_communication(h);
	h.set_storage_type(PST_ADDITIVE);
	m_pFetiLayouts->vec_use_inner_communication(u);
	u.set_storage_type(PST_CONSISTENT);

	// start value
	u.set(0.0);

	// (a) invoke Neumann solver to get \f$h_{\{I \Delta\}}^{(p)}\f$
	if(!m_pNeumannSolver->apply_return_defect(u, h))
	{
		UG_LOG_ALL_PROCS("ERROR in 'PrimalSubassembledMatrixInverse::apply': "
						 "Could not solve Neumann problem (step 2.a) on Proc "
							<< pcl::GetProcRank() << ".\n");
		bSuccess = false;
	}

	write_debug(u, "PSMI_u_2");

	// (b) apply matrix to \f$[h_{\{I \Delta\}}^{(p)}, 0]^T\f$ - multiply with full matrix
	if(!m_pMatrix->apply(h, u))
	{
		UG_LOG_ALL_PROCS("ERROR in 'PrimalSubassembledMatrixInverse::apply': "
						 "Could not apply full matrix (step 2.b) on "
						 "Proc " << pcl::GetProcRank() << ".\n");
		bSuccess = false;
	}

	write_debug(h, "PSMI_h_2_AfterNeumann");

	// (c) compute h = f - h on primal.
	m_pFetiLayouts->vec_scale_add_on_primal(h, 1.0, f, -1.0, h);

	write_debug(h, "PSMI_h_2_Added");

//	2.9 Create storage for u,f on primal root
	vector_type rootF;
	vector_type rootU;
	if(m_primalRootProc == pcl::GetProcRank())
	{
		UG_LOG("Creating coarse vector f of size: " << m_pRootSchurComplementMatrix->num_rows()<<"\n");
		rootF.resize(m_pRootSchurComplementMatrix->num_rows());
		UG_LOG("Creating coarse vector u of size: " << m_pRootSchurComplementMatrix->num_cols()<<"\n");
		rootU.resize(m_pRootSchurComplementMatrix->num_cols());
	}

//	3. Since \f$\tilde{f}_{\Pi}\f$ is saved additively, gather it to one process (root)
//     where it is then consistent.
	rootF.set(0.0);
	VecGather(&rootF, &h, m_masterAllToOneLayout, m_slaveAllToOneLayout);

//	4. Solve \f$S_{\Pi \Pi} u_{\Pi} = \tilde{f}_{\Pi}\f$ on root
//     Toselli, p.~165, below eq.~(6.64): this is a ``local problem with Neumann bnd cnds at edges, zero Dirichlet bnd vrts''

//	only on root proc
	if(m_primalRootProc == pcl::GetProcRank())
	{
	//	only if matrix if non-zero
		if(m_RootSchurComplementOp.get_matrix().num_cols() != 0)
		{
		//	invert matrix
			rootF.set_storage_type(PST_ADDITIVE);
			rootU.set_storage_type(PST_CONSISTENT);
			if(!m_pCoarseProblemSolver->apply(rootU, rootF))
			{
				std::cout << "ERROR in 'PrimalSubassembledMatrixInverse::apply': "
								 "Could not invert Schur complement on root proc."
						<< std::endl;
				bSuccess = false;
			}
		}

	//	Debug output of matrix (only used to flag if debug is wanted)
		if(m_pDebugWriter != NULL)
		{
		//	use own debug writer
			AlgebraDebugWriter<algebra_type, 2> debugWriter;

		//	vector of positions
			std::vector<MathVector<2> > pos(m_pRootSchurComplementMatrix->num_rows());
			if(pos.size()==9)
			{
				pos[0] = MathVector<2>(-0.5, -0.5);
				pos[1] = MathVector<2>(-0.0, -0.5);
				pos[2] = MathVector<2>(-0.5, -0.0);
				pos[3] = MathVector<2>(0.0, 0.0);
				pos[4] = MathVector<2>(0.5, -0.5);
				pos[5] = MathVector<2>(0.5, 0.0);
				pos[6] = MathVector<2>(-0.5, 0.5);
				pos[7] = MathVector<2>(0.0, 0.5);
				pos[8] = MathVector<2>(0.5, 0.5);
			}
		//	set positions
			debugWriter.set_positions(&pos[0], pos.size());

			debugWriter.write_vector(rootU,
										 "RootU");
		}
	}

//	5. Broadcast \f$u_{\Pi}\f$ to all Procs. \f$u_{\Pi}\f$ is consistently saved.
	u.set(0.0);
	VecBroadcast(&u, &rootU, m_slaveAllToOneLayout, m_masterAllToOneLayout);

//	6. Compute other values of solution, keeping values in Primal fixed

	write_debug(u, "PSMI_u_6");

	// (a) set dirichlet values on primal unknowns of rhs
	m_pFetiLayouts->vec_scaled_copy_on_primal(f, u, 1.0);

	write_debug(f, "PSMI_fTmp_6");
 
	// (b) invert neumann problem, with dirichlet values on primal
	m_pFetiLayouts->vec_use_inner_communication(f);
	f.set_storage_type(PST_ADDITIVE);
	m_pFetiLayouts->vec_use_inner_communication(u);
	u.set_storage_type(PST_CONSISTENT);
	if(!m_pNeumannSolver->apply_return_defect(u, f)) // solve with Neumann matrix!
	{
		UG_LOG_ALL_PROCS("ERROR in 'PrimalSubassembledMatrixInverse::apply': "
						 "Could not solve Neumann problem (step 7) on Proc "
							<< pcl::GetProcRank() << ".\n");
		bSuccess = false;
	}

	write_debug(u, "PSMI_u_8");

//	check all procs
	if(!pcl::AllProcsTrue(bSuccess))
	{
		UG_LOG("ERROR in PrimalSubassembledMatrixInverse::apply: Some process could not back solve.\n");
		return false;
	}

//	we're done
	return true;
} /* end 'PrimalSubassembledMatrixInverse::apply_return_defect()' */

template <typename TAlgebra>
bool PrimalSubassembledMatrixInverse<TAlgebra>::
apply(vector_type& x, const vector_type& b)
{
	write_debug(b, "PSMI_apply_b");

//	copy defect
	vector_type d; d.resize(b.size());
	d = b;

	write_debug(d, "PSMI_apply_d");

//	solve on copy of defect
	return apply_return_defect(x, d);
} /* end 'PrimalSubassembledMatrixInverse::apply()' */


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

//	3.1 set layouts in PrimalSubassembledMatrixInverse
	m_PrimalSubassembledMatrixInverse.set_feti_layouts(m_fetiLayouts);

//	3.2 init Neumann system and solver
//	check that neumann solver has been set
	if(m_pNeumannSolver == NULL)
	{
		UG_LOG("ERROR in FETISolver::init: No neumann solver set "
				" for inversion of A_{I,Delta}{I,Delta} in PrimalSubassembledMatrixInverse.\n");
		return false;
	}

//	set neumann solver in PrimalSubassembledMatrixInverse
	m_PrimalSubassembledMatrixInverse.set_neumann_solver(*m_pNeumannSolver);

//  3.3 init coarse problem solver used in PrimalSubassembledMatrixInverse
//	check that coarse problem solver has been set
	if(m_pCoarseProblemSolver == NULL)
	{
		UG_LOG("ERROR in FETISolver::init: No coarse problem solver set "
				" for solving S_{Pi Pi} u_{Pi} = tilde{f}_{Pi}.\n");
		return false;
	}

//	set coarse problem solver in PrimalSubassembledMatrixInverse
	m_PrimalSubassembledMatrixInverse.set_coarse_problem_solver(*m_pCoarseProblemSolver);

//	3.4 init PrimalSubassembledMatrixInverse (operator - given as parameter here - is also set thereby)
	if(m_PrimalSubassembledMatrixInverse.init(*m_pOperator) != true)
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

	bool bSuccess = true;

//	Solve: A u = f
	if(!m_PrimalSubassembledMatrixInverse.apply(u, f))
	{
		UG_LOG("ERROR in FETISolver::apply: Cannot back solve.\n");
		bSuccess = false;
	}

	write_debug(u, "FETI_Sol_Final");

//	check all procs
	if(!pcl::AllProcsTrue(bSuccess))
	{
		UG_LOG("ERROR in FETISolver::apply: Some process could not back solve.\n");
		return false;
	}

	return m_pConvCheck->post();
} /* end 'FETISolver::apply_return_defect()' */

template <typename TAlgebra>
bool FETISolver<TAlgebra>::
apply_F(vector_type& f, const vector_type& v)
{
//	Help vector
	vector_type fTmp; fTmp.create(v.size());
	fTmp.set_storage_type(PST_CONSISTENT);

//	0. Reset values of f
	f.set(0.0);

//	1. Apply transposed jump operator: f = B_{\Delta}^T * v_{\Delta}:
	ComputeDifferenceOnDeltaTransposed(f, v, m_fetiLayouts.get_dual_master_layout(),
	                                   	   	 m_fetiLayouts.get_dual_slave_layout(),
	                                   	   	 m_fetiLayouts.get_dual_nbr_slave_layout());

//  2. Apply PrimalSubassembledMatrixInverse to f
	m_PrimalSubassembledMatrixInverse.apply(fTmp, f);

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

//  1. Apply PrimalSubassembledMatrixInverse to 'f'
	if(!m_PrimalSubassembledMatrixInverse.apply(dTmp, f))
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
template class PrimalSubassembledMatrixInverse<CPUAlgebra>;
template class PrimalSubassembledMatrixInverse<CPUBlockAlgebra<3> >;
template class FETISolver<CPUAlgebra>;
template class FETISolver<CPUBlockAlgebra<3> >;
}// end of namespace



/*	Temporaerer code, der in PrimalSubassembledMatrixInverse eingebaut werden muss.
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
