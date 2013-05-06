// 01.02.2011 (m,d,y)

//	THIS FILE IS ONLY TEMPORARY!
//	I (Sebastian) created this file to reduce compile-time during
//	the creation of the FETI solver. When everything works as
//	expected one should move the code in this file to
//	feti_impl.h, so that it will work with all the different algebras.
//
//	In the moment template instantiations are invoked at the end of this file.

#ifdef UG_PARALLEL

// extern headers
#include <cmath>
#include <sstream>  // added for 'stringstream'

// own header
#include "feti.h"

// algebra types
#include "lib_algebra/cpu_algebra_types.h"
#include "lib_algebra/operator/algebra_debug_writer.h"

// additions for profiling
#include "common/profiler/profiler.h"
#define PROFILE_FETI
#ifdef PROFILE_FETI
	#define FETI_PROFILE_FUNC()			PROFILE_FUNC_GROUP("algebra feti")
	#define FETI_PROFILE_BEGIN(name)	PROFILE_BEGIN_GROUP(name, "algebra feti")
	#define FETI_PROFILE_END()			PROFILE_END()
#else
	#define FETI_PROFILE_FUNC()
	#define FETI_PROFILE_BEGIN(name)
	#define FETI_PROFILE_END()
#endif
// additions for profiling - end

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

template <int dim>
struct PosAndIndex{
		PosAndIndex() {}
		PosAndIndex(int i1, const MathVector<dim>& val) :
		ind1(i1), pos(val) {}
	int ind1;
	MathVector<dim> pos;
};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	LocalSchurComplement implementation
template <typename TAlgebra>
LocalSchurComplement<TAlgebra>::
LocalSchurComplement() :
	m_pMatrix(NULL),
	m_pFetiLayouts(NULL),
	m_spDirichletOperator(new MatrixOperator<matrix_type, vector_type>),
	m_pDirichletMatrix(NULL),
	m_spDirichletSolver(NULL)
{
}

template <typename TAlgebra>
void LocalSchurComplement<TAlgebra>::
init()
{
//	check that operator has been set
	if(m_spOperator.invalid())
		UG_THROW("LocalSchurComplement::init: No Operator A set.");

//	check Feti layouts have been set
	if(m_pFetiLayouts == NULL)
		UG_THROW("LocalSchurComplement::init: FETI layouts not set.");

//	save matrix from which we build the Schur complement
	m_pMatrix = &m_spOperator->get_matrix();

//	get matrix from dirichlet operator
	m_pDirichletMatrix = &m_spDirichletOperator->get_matrix();

//	Copy Matrix for Dirichlet Problem
	*m_pDirichletMatrix = *m_pMatrix;

//	Set Dirichlet values on Pi
	m_pFetiLayouts->mat_set_dirichlet_on_primal(*m_pDirichletMatrix);

//	Set Dirichlet values on Delta
	m_pFetiLayouts->mat_set_dirichlet_on_dual(*m_pDirichletMatrix);

//	Let Dirichlet Matrix use intra subdomain layouts
	m_pFetiLayouts->mat_use_intra_sd_communication(*m_pDirichletMatrix);

//	init sequential solver for Dirichlet problem
	if(m_spDirichletSolver.valid())
		if(!m_spDirichletSolver->init(m_spDirichletOperator))
			UG_THROW("LocalSchurComplement::init: Cannot init "
					"Sequential Dirichlet Solver for Operator A.");

//	Debug output of matrices
	write_debug(m_spDirichletOperator->get_matrix(), "FetiDirichletMatrix.mat");
	write_debug(m_spOperator->get_matrix(), "FetiOriginalMatrix.mat");

//	reset apply counter
	m_applyCnt = 0;

} /* end 'LocalSchurComplement::init()' */

template <typename TAlgebra>
void LocalSchurComplement<TAlgebra>::
apply(vector_type& f, const vector_type& u)
{
	FETI_PROFILE_BEGIN(FETI_apply);
//	check that matrix has been set
	if(m_spOperator.invalid())
		UG_THROW("LocalSchurComplement::apply: Matrix A not set.");

//	check Dirichlet solver
	if(m_spDirichletSolver.invalid())
		UG_THROW("LocalSchurComplement::apply: No sequential Dirichlet Solver set.");

//	Check parallel storage type of matrix
	if(!m_pDirichletMatrix->has_storage_type(PST_ADDITIVE))
		UG_THROW("LocalSchurComplement::apply: Inadequate storage format of matrix.");

//	Check parallel storage type of vectors
	if (!u.has_storage_type(PST_CONSISTENT))
		UG_THROW("LocalSchurComplement::apply: Inadequate storage format of Vector 'u' (should be consistent).");

	if(!f.has_storage_type(PST_ADDITIVE))
		UG_THROW("LocalSchurComplement::apply: Inadequate storage format of Vector 'f' (should be additive).");

//	Help vector
//	\todo: it would be sufficient to copy only the layouts without copying the values
	vector_type uTmp; uTmp.create(u.size()); uTmp = u;

//	1. Set values to zero on I and Pi
	// (a) Reset all values
	uTmp.set(0.0);

	// (b) Copy values on \Delta
	m_pFetiLayouts->vec_scale_append_on_dual(uTmp, u, 1.0);

//	2. Compute rhs f_{I} = A_{I \Delta} u_{\Delta}
//	f is additive afterwards
	m_spDirichletOperator->apply(f, uTmp);

	// set values to zero on \Delta (values are already zero on primal after
	// application of DirichletOperator!)
	m_pFetiLayouts->vec_set_on_dual(f, 0.0);

//	Debug output of vector
	if(debug_writer() != NULL)
	{
	//	add iter count to name
		std::string name("FetiDirichletRhs");
		char ext[20]; sprintf(ext, "_apply%03d.vec", m_applyCnt);
		//name.append(m_statType);
		name.append(ext);

		debug_writer()->write_vector(f, name.c_str());
	}

//	3. Invert on inner unknowns u_{I} = A_{II}^{-1} f_{I}
	// (a) use the intra-FETI-subdomain layouts
	m_pFetiLayouts->vec_use_intra_sd_communication(f);
	f.set_storage_type(PST_ADDITIVE);
	m_pFetiLayouts->vec_use_intra_sd_communication(uTmp);
	uTmp.set_storage_type(PST_CONSISTENT);

	// (b) invoke Dirichlet solver
	//	uTmp is consistent afterwards
	if(!m_spDirichletSolver->apply_return_defect(uTmp, f))
	{
		UG_LOG_ALL_PROCS("ERROR in 'LocalSchurComplement::apply': "
						 "Could not solve Dirichlet problem (step 3.b) on Proc "
							<< pcl::GetProcRank() << " (m_statType = '" << m_statType << "').\n");
		UG_LOG_ALL_PROCS("ERROR in 'LocalSchurComplement::apply':"
						" Last defect was " << m_spDirichletSolver->defect() <<
						" after " << m_spDirichletSolver->step() << " steps.\n");

		UG_THROW("Cannot solve Local Schur Complement.");
	} /* else {
		UG_LOG_ALL_PROCS("'LocalSchurComplement::apply':"
						" Last defect after applying Dirichlet solver (step 3.b) was " << m_pDirichletSolver->defect() <<
						" after " << m_pDirichletSolver->step() << " steps.\n");
	}*/

//	remember for statistic
	StepConv stepConv;
	if(!m_statType.empty())
	{
		stepConv.lastDef3b = m_spDirichletSolver->defect();
		stepConv.numIter3b = m_spDirichletSolver->step();

		m_mvStepConv[m_statType].push_back(stepConv);
	}


//	4. Compute result vector
	// (a) Scale u_{I} by -1
	uTmp *= -1.0;

	// (b) Add u_{\Delta} on \Delta
	m_pFetiLayouts->vec_scale_append_on_dual(uTmp, u, 1.0);

	// (c) Multiply with full matrix
	//	f is additive afterwards
	m_spOperator->apply(f, uTmp);

	// make f consistent (on delta is sufficient)
	f.change_storage_type(PST_CONSISTENT);

//	5. Reset all values for I, \Pi
	uTmp = f;
	f.set(0.0);
	m_pFetiLayouts->vec_scale_append_on_dual(f, uTmp, 1.0);

//	increase apply counter
	m_applyCnt++;
} /* end 'LocalSchurComplement::apply()' */

template <typename TAlgebra>
void LocalSchurComplement<TAlgebra>::
apply_sub(vector_type& f, const vector_type& u)
{
	FETI_PROFILE_BEGIN(FETI_apply_sub);
//	create new rhs
	vector_type d; d.resize(f.size());

//	solve
	apply(d, u);

//	subtract from vector
	f -= d;
}

template <typename TAlgebra>
void LocalSchurComplement<TAlgebra>::
print_statistic_of_inner_solver(bool bPrintOnlyAverages) //const
{
	using namespace std;
//	Process Communicator for CommWorld (MPI_WORLD)
	pcl::ProcessCommunicator ProcCom;

	typename map<string, vector<StepConv> >::const_iterator mapIter = m_mvStepConv.begin();

	for(; mapIter != m_mvStepConv.end(); ++mapIter)
	{
	//	write Type
		std::string type                  = (*mapIter).first;
		const vector<StepConv>& vStepConv = (*mapIter).second;
		double avgLastDef3b = -1.0;

		if (!bPrintOnlyAverages) {
			UG_LOG("Calls of Dirichlet solver in 'LocalSchurComplement::apply' for '"<< type << "' ('avg' is average over procs):\n");
		} else {
			UG_LOG("dps calls for '"<< type << "' (3b): ");
		}

	//	write all calls
	//	print num call header
		if (!bPrintOnlyAverages) {
			UG_LOG("Call                     :  ");
			for(size_t i = 0; i < vStepConv.size(); ++i)
				UG_LOG(std::setw(8) << i << " |  ");
			UG_LOG("\n");
		}

		if (!bPrintOnlyAverages) UG_LOG("Defect3b  (avg)          :  ");

	//	print results
		for(size_t i = 0; i < vStepConv.size(); ++i)
		{
			double tGlob, tLoc = vStepConv[i].lastDef3b;
			ProcCom.allreduce(&tLoc, &tGlob, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
			tGlob /= pcl::GetNumProcesses();
			avgLastDef3b = tGlob;
			if (!bPrintOnlyAverages) UG_LOG(std::setprecision(2) << tGlob << " |  ");
		}
		if (!bPrintOnlyAverages) {
			UG_LOG("\n");

			UG_LOG("NumIter3b (avg, max, min):");
		}
		int sumAvgNumIter = 0;
		for(size_t i = 0; i < vStepConv.size(); ++i)
		{
			int tGlob, tLoc = vStepConv[i].numIter3b;
			ProcCom.allreduce(&tLoc, &tGlob, 1, PCL_DT_INT, PCL_RO_SUM);
			tGlob /= pcl::GetNumProcesses();

			if (!bPrintOnlyAverages) UG_LOG(std::setw(3) << tGlob << ",");

			sumAvgNumIter += tGlob;
			m_totalIterCntOfInnerSolvers += tGlob; // sum over all calls of inner solvers

			tLoc = vStepConv[i].numIter3b;
			ProcCom.allreduce(&tLoc, &tGlob, 1, PCL_DT_INT, PCL_RO_MAX);

			if (!bPrintOnlyAverages) UG_LOG(std::setw(3) << tGlob << ",");

			tLoc = vStepConv[i].numIter3b;
			ProcCom.allreduce(&tLoc, &tGlob, 1, PCL_DT_INT, PCL_RO_MIN);
			if (!bPrintOnlyAverages) UG_LOG(std::setw(3) << tGlob << "|");
		}
		if (bPrintOnlyAverages) {
			UG_LOG(" avg last defect: " << std::setprecision(2) << avgLastDef3b << " in ");
		}
		UG_LOG("" << std::setw(3) << sumAvgNumIter << " avg steps" 
//			      << std::fixed << (double)sumAvgNumIter/vStepConv.size() << " per FETI iteration"
//			      << std::scientific << "\n");
			      << " in " << std::setw(3) << vStepConv.size() << " call(s)\n"); // Note: #calls should be always equal #FETI iterations - but one of them is performed before lambda-loop!

		if (!bPrintOnlyAverages) UG_LOG("\n");
	}
} /* end 'LocalSchurComplement::print_statistic_of_inner_solver()' */

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	PrimalSubassembledMatrixInverse implementation
template <typename TAlgebra>
PrimalSubassembledMatrixInverse<TAlgebra>::
PrimalSubassembledMatrixInverse() :
	m_spOperator(NULL),
	m_pMatrix(NULL),
	m_pFetiLayouts(NULL),
	m_spNeumannOperator(new MatrixOperator<matrix_type, vector_type>),
	m_pNeumannMatrix(NULL),
	m_spNeumannSolver(NULL),
	m_spCoarseProblemSolver(NULL),
	m_primalRootProc(-1),
	m_spRootSchurComplementOp(new MatrixOperator<matrix_type, vector_type>),
	m_pRootSchurComplementMatrix(NULL),
	m_statType(""),
	m_bTestOneToManyLayouts(false)
{
}

template <typename TAlgebra>
bool PrimalSubassembledMatrixInverse<TAlgebra>::
init(SmartPtr<ILinearOperator<vector_type> > L)
{
	FETI_PROFILE_BEGIN(FETI_init);
//	status output
	UG_LOG("     % Initializing 'PrimalSubassembledMatrixInverse': \n");

//	success flag
	bool bSuccess = true;

//	remember operator
	m_spOperator = L.template cast_dynamic<MatrixOperator<matrix_type, vector_type> >();

//	check, that operator is correct
	if(m_spOperator.invalid())
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

//	Check all procs - scheint unnoetig! Wird beim Aufrufer, 'FETISolver::init()', am Ende ebenfalls gemacht!
/*
	if(!pcl::AllProcsTrue(bSuccess))
	{
		UG_LOG_ALL_PROCS("ERROR in 'PrimalSubassembledMatrixInverse::init':"
				" Some processes could not init Schur Complement inverse.\n");
		return false;
	}
*/
//	save matrix from which we build the Schur complement (same matrix as
//  'm_pOperator->get_matrix()' in 'LocalSchurComplement::init()')
	m_pMatrix = &m_spOperator->get_matrix();

//	get matrix from Neumann operator
	m_pNeumannMatrix = &m_spNeumannOperator->get_matrix();

//	Copy Matrix for Neumann Problem
	*m_pNeumannMatrix = *m_pMatrix;

//	Set Dirichlet values on Pi
	m_pFetiLayouts->mat_set_dirichlet_on_primal(*m_pNeumannMatrix);

//	Let Neumann Matrix use intra subdomain layouts
	m_pFetiLayouts->mat_use_intra_sd_communication(*m_pNeumannMatrix);

//	status output
	UG_LOG("     %  - initializing 'NeumannSolver' ... ");

//	init sequential solver for Dirichlet problem
	if(m_spNeumannSolver.valid())
		if(!m_spNeumannSolver->init(m_spNeumannOperator))
		{
			UG_LOG("ERROR in 'PrimalSubassembledMatrixInverse::init': Cannot init "
					"Sequential Neumann Solver for Operator A.\n");
			return false;
		}
	UG_LOG("done.\n");

//	Debug output of matrices
	write_debug(m_spNeumannOperator->get_matrix(), "FetiNeumannMatrix.mat");

//	Choose root process, where Schur complement w.r.t. Primal unknowns
//	is gathered.
	m_primalRootProc = 0;

////////////////////////////////////////////////////////////////////////////////
// Nomenclature:
//  'total'   - regarding the whole domain
//  'subdomain' or similar - regarding the feti subdomain the current proc belongs to
//  'local'   - regarding the current proc only
//  'rootIDs' - (local) algebra indices of primal variables on root proc.
////////////////////////////////////////////////////////////////////////////////

//	vector to store newly created root ids - "look-up-table".
//	please note that rootIDs may only be indexed with the algebra-index
//	of primal nodes. Content is not defined for other entries.
	std::vector<int> vPrimalRootIDLUT; // original name: 'rootIDs'

//	Build layouts such that all processes can communicate their unknowns
//	to the primal Root Process
	UG_LOG("     %  - building 'OneToManyLayout' ... ");
	FETI_PROFILE_BEGIN(PrimalSubassMatInvInit_BuildOneToManyLayout);
	int newVecSize = BuildOneToManyLayout(m_masterAllToOneLayout,
						 m_slaveAllToOneLayout, m_primalRootProc,
						 m_pFetiLayouts->get_primal_master_layout(),
						 m_pFetiLayouts->get_primal_slave_layout(),
						 m_allToOneProcessComm, &vPrimalRootIDLUT);
	FETI_PROFILE_END();			// end 'FETI_PROFILE_BEGIN(PrimalSubassMatInvInit_BuildOneToManyLayout)' - Messpunkt ok
	UG_LOG("done.\n");

	if (m_bTestOneToManyLayouts == true) {
		UG_LOG("     %  - TEST ONE TO MANY LAYOUTS:\n");
		pcl::InterfaceCommunicator<IndexLayout> comTmp;
		if (TestLayout(m_spOperator->layouts()->proc_comm(),
				comTmp, m_masterAllToOneLayout, m_slaveAllToOneLayout, true) != true) {
			UG_LOG("     %  - ONE TO MANY LAYOUTS inconsistent!\n");
		} else {
			UG_LOG("     %  - ONE TO MANY LAYOUTS are consistent!\n");
		}
	}

//	We have to gather the primal root IDs and the quantities of primal nodes
//	on each process of the feti subdomain in one array.
//	first we collect all local primal variables in one array
//	Collect all Primal indices on proc
	std::vector<IndexLayout::Element> vLocalPrimalLocalID; // former name: 'vLocalPrimalIndex'
	CollectUniqueElements(vLocalPrimalLocalID, m_slaveAllToOneLayout);

//	now build an array of local primal root ids
	std::vector<int> vLocalPrimalRootID(vLocalPrimalLocalID.size());
	for(size_t i = 0; i < vLocalPrimalLocalID.size(); ++i)
		vLocalPrimalRootID[i] = vPrimalRootIDLUT[vLocalPrimalLocalID[i]];

	pcl::ProcessCommunicator& intraFetiSubdomComm =
						m_pFetiLayouts->get_intra_sd_process_communicator();

//	vector that holds number of primal variables on each process
//	of the local feti subdomain.
	std::vector<int>	vNumPrimalVariablesPerProc;

//	rootIDs of primal variables in the local feti subdomain
	std::vector<int> vSubdomPrimalRootID;

	intraFetiSubdomComm.allgatherv(vSubdomPrimalRootID, vLocalPrimalRootID,
								   &vNumPrimalVariablesPerProc);


//	log num primal variables (display width for moderate sizes of indices;
//  with 'UG_LOG()' only stuff concerning the outproc ...):
	if (intraFetiSubdomComm.size() != vNumPrimalVariablesPerProc.size())
		UG_LOG("     %  - Huh!??: intraFetiSubdomComm.size() = " << intraFetiSubdomComm.size() <<
			   " != vNumPrimalVariablesPerProc.size() = " << vNumPrimalVariablesPerProc.size() << std::endl);

	UG_LOG("     %  -------------------------------------------------------------------" << std::endl); 
	UG_LOG("     %  - Assemble entries of Schur complement matrix on 'primal root proc' " << std::setw(4) << m_primalRootProc   << std::endl);
	UG_LOG("     %  - Log number of primal variables ('#pv') for proc ('output proc')   " << std::setw(4) << pcl::GetProcRank() << std::endl);
	UG_LOG("     %  -    number of processes on this feti subdomain:      " << intraFetiSubdomComm.size() << ","  << std::endl);
	UG_LOG("     %  -                         as looped ([counter:rank]): ");
	for (size_t procInFetiSD = 0; procInFetiSD < intraFetiSubdomComm.size(); ++procInFetiSD)
		UG_LOG("[" << procInFetiSD << ":" << intraFetiSubdomComm.get_proc_id(procInFetiSD) << "] ");
	UG_LOG(std::endl);
	//UG_LOG("     %  -    num procs in 'm_allToOneProcessComm': " << m_allToOneProcessComm.size() << std::endl); // is "world"!


	UG_LOG("     %  -    num primal variables in total (whole domain)");
	if (pcl::GetProcRank() == m_primalRootProc) {
		UG_LOG(": " << newVecSize << " (= 'newVecSize')");
		//for (size_t i = 0; i < vPrimalRootIDLUT.size(); ++i)
		//	UG_LOG(std::setw(6) << vPrimalRootIDLUT[i]);
		UG_LOG(std::endl);
	} else {
		UG_LOG(" only known by primal root proc (try '-outproc <primal root proc rank>')!" << std::endl);
	}

	UG_LOG("     %  -                         on this feti subdomain:  " << vSubdomPrimalRootID.size() << std::endl);
	UG_LOG("     %  -                         on each proc of this subdomain ([rank:#pv]): ");
	for (size_t i = 0; i < vNumPrimalVariablesPerProc.size(); ++i) // 'vNumPrimalVariablesPerProc.size()' = 'intraFetiSubdomComm.size()'!
		UG_LOG("[" << intraFetiSubdomComm.get_proc_id(i) << ":"<< vNumPrimalVariablesPerProc[i] << "] ");
	UG_LOG(std::endl);

	UG_LOG("     %  -    primal root proc algebra indices on this feti subdomain ('vSubdomPrimalRootID's): ");
	for (size_t i = 0; i < vSubdomPrimalRootID.size(); ++i)
		UG_LOG(std::setw(6) << vSubdomPrimalRootID[i] << " ");
	UG_LOG(std::endl);

	UG_LOG("     %  -                                     on proc " <<
		   std::setw(4) << pcl::GetProcRank()<<" (outproc;   'vLocalPrimalRootID's): ");
	for (size_t i = 0; i < vLocalPrimalLocalID.size(); ++i)
		UG_LOG(std::setw(6) << vLocalPrimalRootID[i] << " "); // = 'vPrimalRootIDLUT[vLocalPrimalLocalID[i]]' (see above)
	UG_LOG(std::endl);

	UG_LOG("     %  -                   corresponding local-proc algebra indices ('vLocalPrimalLocalID's): ");
	for (size_t i = 0; i < vLocalPrimalLocalID.size(); ++i)
		UG_LOG(std::setw(6) << vLocalPrimalLocalID[i] << " ");
	UG_LOG(std::endl);

// Some checks
	for(size_t procInFetiSD = 0; procInFetiSD < intraFetiSubdomComm.size(); procInFetiSD++) {
		if(pcl::GetProcRank() == intraFetiSubdomComm.get_proc_id(procInFetiSD))
			if (vNumPrimalVariablesPerProc[procInFetiSD] != (int)vLocalPrimalLocalID.size()) {
				UG_LOG("     %  - Huh? 'vNumPrimalVariablesPerProc[" << pcl::GetProcRank() << "]' != 'vLocalPrimalLocalID.size()': ");
				UG_LOG(vNumPrimalVariablesPerProc[procInFetiSD] << " != " << vLocalPrimalLocalID.size() << "!?\n");
			} //else {
		//UG_LOG("     %  -    'vNumPrimalVariablesPerProc[" << pcl::GetProcRank() << "]' == 'vLocalPrimalLocalID.size()': ");
		//UG_LOG(vNumPrimalVariablesPerProc[procInFetiSD] << " == " << vLocalPrimalLocalID.size() << "\n");
		//}
	}

//	processes will collect their local primal connections here.
	typedef PrimalConnection<typename vector_type::value_type> PrimalConnection;
	std::vector<PrimalConnection> vLocalPrimalConnections;

//	create help vectors
	vector_type e;   e.resize(m_pMatrix->num_rows());
	vector_type h1; h1.resize(m_pMatrix->num_rows());
	vector_type h2; h2.resize(m_pMatrix->num_rows());

//	set communication to intra subdomain communication
	m_pFetiLayouts->vec_use_intra_sd_communication(e);
	m_pFetiLayouts->vec_use_intra_sd_communication(h1);
	m_pFetiLayouts->vec_use_intra_sd_communication(h2);

	UG_LOG("     %\n");
	UG_LOG("     %  - assemble entries of 'S_PiPi' (format: '[<proc rank>]: (<from_i> -> <to_j>)') ... \n");

	FETI_PROFILE_BEGIN(PrimalSubassMatInvInit_Assemble_S_PiPi);
//	Now within each feti subdomain, the primal unknowns are looped one after the
//	other. This is done like the following: Each process loops the number of
//	procs of the feti subdomain it belongs, and then over their primal unknowns
//	(as stored in 'vNumPrimalVariablesPerProc') on this subdomain. That way
//	for every primal unknown on the feti subdomain, we can set its value to 1
//	while all the other values are zero. This gives us the unity vector. For
//	this unity vector an application of S_{Pi Pi} is computed such that we can
//	read then the values of (S_{Pi Pi})_{ij} in the i-th component of the result
//	vector. All those couplings are stored in the vector of connections called
//	vLocalPrimalConnections and sent to the primalRootProc at the end of the
//	loop.
	size_t primalCounter = 0;
	std::stringstream ss;

	for(size_t procInFetiSD = 0; procInFetiSD < intraFetiSubdomComm.size();
			procInFetiSD++)
	{
		UG_LOG("     %  - [proc " << std::setw(6) << intraFetiSubdomComm.get_proc_id(procInFetiSD) << "]: ");

		// loop over feti subdomain primals (\Pi variables) of current proc in loop
		for(size_t pvTo_j = 0; pvTo_j < (size_t)vNumPrimalVariablesPerProc[procInFetiSD]; ++pvTo_j)
		{
		////////////////////////////////////////
		// 	1. Create (column) unity vector e_j:
		////////////////////////////////////////

		//	keep the matrix id on root, of the primal unknown, that is set to one
			int unityRootID = vSubdomPrimalRootID[primalCounter++];

		//	reset identity vector to zero for all (primal) unknowns
			e.set(0.0);

			ss.str(""); // clear stringstream before next loop

		//	set value of unity vector to one if on process and quantity, else 0
			if(pcl::GetProcRank() == intraFetiSubdomComm.get_proc_id(procInFetiSD))
			{
				const IndexLayout::Element localPrimalIndex = vLocalPrimalLocalID[pvTo_j];

				e[localPrimalIndex] = 1.0;
				ss << std::setw(3) << vPrimalRootIDLUT[localPrimalIndex]; // store for later output
			} else {
				ss << "no connection to pvTo_j " << pvTo_j;
				//continue; // as long as null space of A_{\{I, \Delta\} \Pi} is empty
			}

		////////////////////////////////////////////////////////////////////////
		// 	2. Apply first matrix A_{\{I, \Delta\} \Pi} to unity vector e^{(p)}:
		////////////////////////////////////////////////////////////////////////
		//	build h1 = A_{\{I \Delta\} \Pi} e
			m_spOperator->apply(h1, e);

		///////////////////////////////////////////////////////////////////
		//	3. Apply A_{\{I \Delta\}\{I\Delta\}}^{-1} by solving "I,\Delta"
		//	   subsystem problem:
		///////////////////////////////////////////////////////////////////
		//	Solve: A_{\{I \Delta\}  \{I \Delta\} } h2 = h1

		//	(a1) Set zero dirichlet bnd conds for rhs h1 on_\Pi
			m_pFetiLayouts->vec_set_on_primal(h1, 0.0);

		//	(a2) Start with zero iterate (not obligatory)
			h2.set(0.0);

		//	(b) Solve dirichlet problem. (Dirichlet rows in A for \Pi dof's are
		//	    already set in 'SchurComplementInverse::init()'!)
			m_pFetiLayouts->vec_use_intra_sd_communication(h1);
			h1.set_storage_type(PST_ADDITIVE);
			m_pFetiLayouts->vec_use_intra_sd_communication(h2);
			h2.set_storage_type(PST_CONSISTENT);

			FETI_PROFILE_BEGIN(PSMIInit_NeumannSolve_SC);
			if(!m_spNeumannSolver->apply(h2, h1))
			{
				UG_LOG_ALL_PROCS("ERROR in 'PrimalSubassembledMatrixInverse::init':"
						" Could not solve local Neumann problem (inversion of A_I Delta,I Delta)"
						" to compute Schur complement w.r.t. primal unknowns: "
								 << procInFetiSD << ", " << pvTo_j << ".\n"); // TODO: info about current repetition of loop? (28102011ih)

				UG_LOG_ALL_PROCS("ERROR in 'PrimalSubassembledMatrixInverse::init':"
								" Last defect was " << m_spNeumannSolver->defect() <<
								" after " << m_spNeumannSolver->step() << " steps.\n");

				return false;
			}
			FETI_PROFILE_END(); // end 'FETI_PROFILE_BEGIN(PSMIInit_NeumannSolve_SC)'

		//	remember for statistic
			StepConv stepConv;
			if(!m_statType.empty())
			{
				stepConv.lastDefSC = m_spNeumannSolver->defect();
				stepConv.numIterSC = m_spNeumannSolver->step();

				m_mvStepConv[m_statType].push_back(stepConv);
			}

		//////////////////////////
		// 	4. Apply third matrix A_{\Pi \{I, \Delta\}}$ to $h_2^{(p)}
		//////////////////////////

		//	(a) Set h2 zero on \Pi. This is enforced by neumann solver

		//	(b) Apply third matrix: h1 = A h2
			m_spOperator->apply(h1, h2);

		//	(c) Set entries to zero on I, \Delta (not needed, therefore skipped)

		///////////////////////////
		// 	5. Compute first term, application of A_{\Pi \Pi} on unity vector
		///////////////////////////

		//	(a) multiply unity vector with matrix h2 = A e
			m_spOperator->apply(h2, e);

		//	(b) Set entries to zero on I, \Delta (not needed, therefore skipped)

		///////////////////////////
		// 	6. Add parts to get result
		///////////////////////////

		//	e = h2 - h1
			m_pFetiLayouts->vec_scale_add_on_primal(e, 1.0, h2, -1.0, h1);


		// 	at this point, we have the contribution of S_ij^{p} in all primal
		//	unknowns i. Thus, we have to read it and send it to the root process

		//	loop process local primal unknowns
			for(size_t pvFrom_i = 0; pvFrom_i < vLocalPrimalLocalID.size(); ++pvFrom_i)
			{
				const IndexLayout::Element localPrimalIndex = vLocalPrimalLocalID[pvFrom_i];

				UG_LOG("(" << std::setw(3) << vPrimalRootIDLUT[localPrimalIndex] <<
					   " -> " << ss.str() << ") ");

			//	read coupling
				typename vector_type::value_type& entry = e[localPrimalIndex];

			//	read root index
				int primalRootID = vPrimalRootIDLUT[localPrimalIndex];

			//  remember coupling
				vLocalPrimalConnections.push_back(PrimalConnection(primalRootID,
				                                                  unityRootID, entry));
			}
			if (((size_t)vNumPrimalVariablesPerProc[procInFetiSD] !=0) && vLocalPrimalLocalID.size() != 0)
				UG_LOG(((pvTo_j < (size_t)vNumPrimalVariablesPerProc[procInFetiSD]-1) ?
						"\n     %                   " : "\n"));

		} // end loop over feti subdomain primals 'pvTo_j' of current proc
		if (((size_t)vNumPrimalVariablesPerProc[procInFetiSD] == 0) || vLocalPrimalLocalID.size() == 0)
			UG_LOG((((size_t)vNumPrimalVariablesPerProc[procInFetiSD] != 0) && (vLocalPrimalLocalID.size() != 0) ? "     %+\n" : "\n"));
	} // end loop over procs in feti subdomain

	UG_LOG("     %  - done.\n");
	UG_LOG("     %  -------------------------------------------------------------------" << std::endl); 

// Further checks
	if (primalCounter != vSubdomPrimalRootID.size()) {
		UG_LOG("     %  - Huh? 'primalCounter' != 'vSubdomPrimalRootID.size()': ");
		UG_LOG(primalCounter << " != " << vSubdomPrimalRootID.size() << "!?\n");
	} //else {
		//UG_LOG("     %  -    'primalCounter' == 'vSubdomPrimalRootID.size()': ");
		//UG_LOG(primalCounter << " == " << vSubdomPrimalRootID.size() << "\n");
	//}


//	all processes send their connections to root
// \todo: This could be improved, so that only processes which contain
//		a primal node are involved. - How? 'if (vLocalPrimalLocalID,size() > 0)'??
	std::vector<PrimalConnection> vPrimalConnections;//	only filled on root
	pcl::ProcessCommunicator commWorld;
	commWorld.gatherv(vPrimalConnections, vLocalPrimalConnections, m_primalRootProc);
	FETI_PROFILE_END();			// end 'FETI_PROFILE_BEGIN(PrimalSubassMatInvInit_Assemble_S_PiPi)' - Messpunkt ok

//	build matrix on primalRoot
	if(pcl::GetProcRank() == m_primalRootProc)
	{
		UG_LOG("     %  - On primal root proc: building Schur complement matrix on primal root proc.\n");

	//	get matrix
		m_pRootSchurComplementMatrix = &m_spRootSchurComplementOp->get_matrix();

	//	check matrix
		if(m_pRootSchurComplementMatrix == NULL)
		{
			UG_LOG("ERROR in 'PrimalSubassembledMatrixInverse::init': No matrix in"
					"Root Schur Complement Operator.\n");
			return false;
		}

	//	reference for convenience
		matrix_type& mat = *m_pRootSchurComplementMatrix;

	//	create matrix of correct size
		mat.resize_and_clear(newVecSize, newVecSize);

	//	info output
		UG_LOG("     %  - Creating Schur Complement matrix"
			   " of size '" << newVecSize <<"x"<<newVecSize << "'" <<std::endl);

	//	copy received values into matrix
		mat.set(0.0);
		UG_LOG("     %  - " << vPrimalConnections.size() << " primal connections received" << std::endl);
//		std::cout << "Writing " << vPrimalConnections.size() << " primal connections:" << std::endl;
		for(size_t i = 0; i < vPrimalConnections.size(); ++i)
		{
		//	get sent connection
			PrimalConnection& pc = vPrimalConnections[i];
//			std::cout << "  ind1: " << pc.ind1 << "    ind2: " << pc.ind2 << "    value: " << pc.value << std::endl;

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
		{
			if(m_spCoarseProblemSolver.valid())
			{
				if(!m_spCoarseProblemSolver->init(m_spRootSchurComplementOp))
				{
					UG_LOG("ERROR in 'PrimalSubassembledMatrixInverse::init': Cannot init "
							"coarse problem Solver for Operator S_{Pi Pi}.\n");
					return false;
				}
			}
			else
			{
				UG_LOG("ERROR in 'PrimalSubassembledMatrixInverse::init': S_{Pi Pi} "
						" needs to be inverted, but no CoarseSolver given.\n");
				return false;
			}
		}

	//	set correct parallel storage type of coarse problem matrix (obviously,
	//	this problem is solved in serial by a single process, but the parallel
	//	operations, e.g. computation of defect, require this
		m_pRootSchurComplementMatrix->set_storage_type(PST_ADDITIVE);

	} // end 'if(pcl::GetProcRank() == m_primalRootProc)'

//	Debug output of matrix
//	this is 2d only debug output. \todo: generalize (i.e. copy+paste)
	if(debug_writer() != NULL && debug_writer()->current_dimension() == 2)
	{
	//	vector of root index + pos
		std::vector<PosAndIndex<2> > vProcLocPos;

	//	get positions on local proc
		for(size_t pqi = 0; pqi < vLocalPrimalLocalID.size(); ++pqi)
		{
		//	get local index of primal
			const IndexLayout::Element localPrimalIndex = vLocalPrimalLocalID[pqi];

		//	get position of local index
			MathVector<2> vPos =
					(const_cast<const IDebugWriter<TAlgebra>*>(&(*debug_writer()))->
							template get_positions<2>())[localPrimalIndex];

		//	read root index
			int id = vPrimalRootIDLUT[localPrimalIndex];

		//	add position to send buffer
			vProcLocPos.push_back(PosAndIndex<2>(id, vPos));
		}

	//	gather positions to root
		std::vector<PosAndIndex<2> > vPosRootSchur;//	only filled on root
		pcl::ProcessCommunicator commWorld;
		commWorld.gatherv(vPosRootSchur, vProcLocPos, m_primalRootProc);

	//	only root proc writes the matrix
		if(pcl::GetProcRank() == m_primalRootProc)
		{
		//	positions are now stored as in vPrimalRootIDLUT, sort this
			std::vector<MathVector<2> > vPosRootSchurSorted(newVecSize);
			for(size_t i = 0; i < vPosRootSchur.size(); ++i)
				vPosRootSchurSorted[vPosRootSchur[i].ind1] = vPosRootSchur[i].pos;

		//	create algebra debug writer
			AlgebraDebugWriter<algebra_type, 2> dbgWriter2d;
			dbgWriter2d.set_positions(&vPosRootSchurSorted[0], newVecSize);

		//	write matrix
			dbgWriter2d.write_matrix(m_spRootSchurComplementOp->get_matrix(),
										 "RootSchurComplementMatrix.mat");
		}
	}


	(void) bSuccess; // removes unused variable-warning
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
	if(m_spNeumannSolver.invalid())
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
	// use intra subdomain interfaces for solving
	m_pFetiLayouts->vec_use_intra_sd_communication(h);
	h.set_storage_type(PST_ADDITIVE);
	m_pFetiLayouts->vec_use_intra_sd_communication(u);
	u.set_storage_type(PST_CONSISTENT);

	// start value
	u.set(0.0);

	// (a) invoke Neumann solver to get \f$u_{\{I \Delta\}}^{(p)}\f$
	FETI_PROFILE_BEGIN(PSMIApply_NeumannSolve_2a);
	if(!m_spNeumannSolver->apply_return_defect(u, h))
	{
		UG_LOG_ALL_PROCS("ERROR in 'PrimalSubassembledMatrixInverse::apply': "
						 "Could not solve Neumann problem (step 2.a) on Proc "
							<< pcl::GetProcRank() << " (m_statType = '" << m_statType << "').\n");

		UG_LOG_ALL_PROCS("ERROR in 'PrimalSubassembledMatrixInverse::apply':"
						" Last defect was " << m_spNeumannSolver->defect() <<
						" after " << m_spNeumannSolver->step() << " steps.\n");
		bSuccess = false;
	}
	FETI_PROFILE_END(); // end 'FETI_PROFILE_BEGIN(PSMIApply_NeumannSolve_2a)'

//	remember for statistic
	StepConv stepConv;
	if(!m_statType.empty())
	{
		stepConv.lastDef2a = m_spNeumannSolver->defect();
		stepConv.numIter2a = m_spNeumannSolver->step();
	}

//	save current solution - 'u' is overwritten by broadcasting \f$u_{\Pi}\f$ after solving (21022011ih)
	vector_type uTmp;  uTmp.create(u.size());
	uTmp = u;

	// (b) apply matrix to \f$[u_{\{I \Delta\}}^{(p)}, 0]^T\f$ - multiply with full matrix
	if(!m_pMatrix->apply(h, u))
	{
		UG_LOG_ALL_PROCS("ERROR in 'PrimalSubassembledMatrixInverse::apply': "
						 "Could not apply full matrix (step 2.b) on "
						 "Proc " << pcl::GetProcRank() << ".\n");
		bSuccess = false;
	}

	// (c) compute h = f - h on primal (this h corresponds to \f$\tilde{f}_{\Pi}^{(p)}\f$!)
	m_pFetiLayouts->vec_scale_add_on_primal(h, 1.0, f, -1.0, h);

//	Create storage for u,f on primal root
	vector_type rootF;
	vector_type rootU;
	if(m_primalRootProc == pcl::GetProcRank())
	{
//		UG_LOG("Creating coarse vector f of size: " << m_pRootSchurComplementMatrix->num_rows()<<"\n");
		rootF.resize(m_pRootSchurComplementMatrix->num_rows());
//		UG_LOG("Creating coarse vector u of size: " << m_pRootSchurComplementMatrix->num_cols()<<"\n");
		rootU.resize(m_pRootSchurComplementMatrix->num_cols());
	}

//	3. Since \f$\tilde{f}_{\Pi}\f$ is saved additively, gather it to one process (root)
//     where it is then consistent.
	rootF.set(0.0);
	//pcl::SynchronizeProcesses();			// TMP
	FETI_PROFILE_BEGIN(PSMIApply_VecGather);
	VecGather(&rootF, &h, m_masterAllToOneLayout, m_slaveAllToOneLayout);
	FETI_PROFILE_END(); // end 'FETI_PROFILE_BEGIN(PSMIApply_VecGather)'

//	4. Solve \f$S_{\Pi \Pi} u_{\Pi} = \tilde{f}_{\Pi}\f$ on root
//     Toselli, p.~165, below eq.~(6.64): this is a ``local problem with Neumann bnd cnds at edges, zero Dirichlet bnd vrts''

//	only on root proc
	if(m_primalRootProc == pcl::GetProcRank())
	{
	//	only if matrix is non-zero
		if(m_spRootSchurComplementOp->get_matrix().num_cols() != 0)
		{
		//	invert matrix
			rootF.set_storage_type(PST_ADDITIVE);
			rootU.set_storage_type(PST_CONSISTENT);

			FETI_PROFILE_BEGIN(PSMIApply_SolveCoarseProblem);
			if(!m_spCoarseProblemSolver->apply_return_defect(rootU, rootF))
			{
				std::cout << "ERROR in 'PrimalSubassembledMatrixInverse::apply': "
								 "Could not invert Schur complement 'S_PiPi' on root proc."
						<< std::endl;
				bSuccess = false;
			} /*
			else {			// TMP
					UG_LOG("'PrimalSubassembledMatrixInverse::apply':"
						   " Last defect after applying coarse problem solver (step 4  ) was " << m_pCoarseProblemSolver->defect() <<
						   " after " << m_pCoarseProblemSolver->step() << " steps.\n");
			}*/

			FETI_PROFILE_END();	// end 'FETI_PROFILE_BEGIN(PSMIApply_SolveCoarseProblem)' - Messpunkt ok, da nur auf einem Proc
		}
	}

//	5. Broadcast \f$u_{\Pi}\f$ to all Procs. \f$u_{\Pi}\f$ is consistently saved.
	u.set(0.0);
	VecBroadcast(&u, &rootU, m_slaveAllToOneLayout, m_masterAllToOneLayout);

//	6.  create help vectors
	vector_type t;  t.create(u.size());
	vector_type uPi; uPi.create(u.size());

	// (a) Copy \f$u_{\Pi}\f$
	uPi.set(0.0);
	m_pFetiLayouts->vec_scale_assign_on_primal(uPi, u, 1.0);

	// (b) apply matrix to \f$[0, u_{\Pi}^{(p)}]^T\f$ - multiply with full matrix
	if(!m_pMatrix->apply(t, uPi))
	{
		UG_LOG_ALL_PROCS("ERROR in 'PrimalSubassembledMatrixInverse::apply': "
						 "Could not apply full matrix (step 5.1) on "
						 "Proc " << pcl::GetProcRank() << ".\n");
		bSuccess = false;
	}

	// (c) Reset values on Pi
	m_pFetiLayouts->vec_set_on_primal(t, 0.0);

//	7. Solve (again) \f$A_D \cdot u = t\f$, using start value \f$u = 0\f$
	vector_type uTmp2;  uTmp2.create(u.size());
	uTmp2.set(0.0);

	// (b) invert neumann problem, with dirichlet values on primal
	m_pFetiLayouts->vec_use_intra_sd_communication(t);
	t.set_storage_type(PST_ADDITIVE);
	m_pFetiLayouts->vec_use_intra_sd_communication(uTmp2);
	uTmp2.set_storage_type(PST_CONSISTENT);

	FETI_PROFILE_BEGIN(PSMIApply_NeumannSolve_7);
	if(!m_spNeumannSolver->apply_return_defect(uTmp2, t)) // solve with Neumann matrix!
	{
		UG_LOG_ALL_PROCS("ERROR in 'PrimalSubassembledMatrixInverse::apply': "
						 "Could not solve Neumann problem (step 7) on Proc "
							<< pcl::GetProcRank() << " (m_statType = '" << m_statType << "').\n");

		UG_LOG_ALL_PROCS("ERROR in 'PrimalSubassembledMatrixInverse::apply':"
						" Last defect was " << m_spNeumannSolver->defect() <<
						" after " << m_spNeumannSolver->step() << " steps.\n");
		bSuccess = false;
	}
	FETI_PROFILE_END(); // end 'FETI_PROFILE_BEGIN(PSMIApply_NeumannSolve_7)'

//	remember for statistic
	if(!m_statType.empty())
	{
		stepConv.lastDef7 = m_spNeumannSolver->defect();
		stepConv.numIter7 = m_spNeumannSolver->step();

		m_mvStepConv[m_statType].push_back(stepConv);
	}

	// (c) compute u = uTmp - uTmp2 on inner and dual
	m_pFetiLayouts->vec_set_on_primal(uTmp2, 0.0);
	u = uTmp;
	u -= uTmp2;

//	assemble solution for primal variables (already computed earlier) into complete solution
	m_pFetiLayouts->vec_scale_assign_on_primal(u, uPi, 1.0);

//	check all procs
	if(!pcl::AllProcsTrue(bSuccess))
	{
		UG_LOG("ERROR in PrimalSubassembledMatrixInverse::apply: Some processes "
			   "could not apply 'PrimalSubassembledMatrixInverse' (end of method).\n");
		return false;
	}

//	we're done
	return true;
} /* end 'PrimalSubassembledMatrixInverse::apply_return_defect()' */

template <typename TAlgebra>
bool PrimalSubassembledMatrixInverse<TAlgebra>::
apply(vector_type& x, const vector_type& b)
{
	FETI_PROFILE_FUNC();
//	copy defect
	vector_type d; d.resize(b.size());
	d = b;

//	solve on copy of defect
	return apply_return_defect(x, d);
} /* end 'PrimalSubassembledMatrixInverse::apply()' */


template <typename TAlgebra>
void PrimalSubassembledMatrixInverse<TAlgebra>::
print_statistic_of_inner_solver(bool bPrintOnlyAverages) //const
{
	using namespace std;
//	Process Communicator for CommWorld (MPI_WORLD)
	pcl::ProcessCommunicator ProcCom;

	typename map<string, vector<StepConv> >::const_iterator mapIter = m_mvStepConv.begin();

	if (!bPrintOnlyAverages) {
		UG_LOG("Calls of Neumann solver in 'PrimalSubassembledMatrixInverse::init' and 'apply_return_defect' ('avg' is average over procs):\n");
	}

	for(; mapIter != m_mvStepConv.end(); ++mapIter)
	{
	//	write Type
		std::string type                  = (*mapIter).first;
		const vector<StepConv>& vStepConv = (*mapIter).second;
		double avgLastDefSC, avgLastDef2a, avgLastDef7;
		avgLastDefSC = avgLastDef2a = avgLastDef7 = -1.0;
		int sumAvgNumIter = -1;

	//	write all calls
	//	print results
		if (type.compare("PSMI_init") == 0) {

			if (!bPrintOnlyAverages) {
				UG_LOG("Calls of Neumann solver in 'PrimalSubassembledMatrixInverse::init' for '"<< type << "':\n");
			} else {
				UG_LOG("nps calls for '"<< type << "'         (SC): ");
			}

	//	print num call header
			if (!bPrintOnlyAverages) {
				UG_LOG("Call                     :  ");
				for(size_t i = 0; i < vStepConv.size(); ++i)
					UG_LOG(std::setw(8) << i << " |  ");
				UG_LOG("\n");
			}

			if (!bPrintOnlyAverages) UG_LOG("DefectSC  (avg)          :  ");

	// (no iteration for this tag!)

			double tGlob, tLoc = vStepConv[0].lastDefSC;
			ProcCom.allreduce(&tLoc, &tGlob, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
			tGlob /= pcl::GetNumProcesses();
			avgLastDefSC = tGlob;
			if (!bPrintOnlyAverages) {
				UG_LOG(std::setprecision(2) << tGlob << " |  ");
				UG_LOG("\n");

				UG_LOG("NumIterSC (avg, max, min):");
			}

			sumAvgNumIter = 0;
			int tGlobN, tLocN = vStepConv[0].numIterSC;
			ProcCom.allreduce(&tLocN, &tGlobN, 1, PCL_DT_INT, PCL_RO_SUM);
			tGlobN /= pcl::GetNumProcesses();

			if (!bPrintOnlyAverages) UG_LOG(std::setw(3) << tGlobN << ",");

			sumAvgNumIter += tGlobN;
			m_totalIterCntOfInnerSolvers += tGlobN; // sum over all calls of inner solvers

			tLocN = vStepConv[0].numIterSC;
			ProcCom.allreduce(&tLocN, &tGlobN, 1, PCL_DT_INT, PCL_RO_MAX);

			if (!bPrintOnlyAverages) UG_LOG(std::setw(3) << tGlobN << ",");

			tLocN = vStepConv[0].numIterSC;
			ProcCom.allreduce(&tLocN, &tGlobN, 1, PCL_DT_INT, PCL_RO_MIN);

			if (!bPrintOnlyAverages) UG_LOG(std::setw(3) << tGlobN << "|");

			if (bPrintOnlyAverages) {
				UG_LOG(" avg last defect: " << std::setprecision(2) << avgLastDefSC << " in ");
			}
			UG_LOG("" << std::setw(3) << sumAvgNumIter << " avg steps\n");

		} else {

		if (!bPrintOnlyAverages) {
			UG_LOG("Calls of Neumann solver in 'PrimalSubassembledMatrixInverse::apply_return_defect' for '"<< type << "'");
			if (type.compare("compute_d") == 0) UG_LOG(" (1x before lambda loop)");
			if (type.compare("backsolve") == 0) UG_LOG(" (1x after  lambda loop)");
			UG_LOG(":\n");
		} else {
			UG_LOG("nps calls for '"<< type << "'         (2a): ");
		}

	//	print num call header
		if (!bPrintOnlyAverages) {
			UG_LOG("Call                     :  ");
			for(size_t i = 0; i < vStepConv.size(); ++i)
				UG_LOG(std::setw(8) << i << " |  ");
			UG_LOG("\n");
		}

		if (!bPrintOnlyAverages) UG_LOG("Defect2a  (avg)          :  ");

	//	print results
		for(size_t i = 0; i < vStepConv.size(); ++i)
		{
			double tGlob, tLoc = vStepConv[i].lastDef2a;
			ProcCom.allreduce(&tLoc, &tGlob, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
			tGlob /= pcl::GetNumProcesses();
			avgLastDef2a = tGlob;
			if (!bPrintOnlyAverages) UG_LOG(std::setprecision(2) << tGlob << " |  ");
		}
		if (!bPrintOnlyAverages) {
			UG_LOG("\n");

			UG_LOG("NumIter2a (avg, max, min):");
		}
		sumAvgNumIter = 0;
		for(size_t i = 0; i < vStepConv.size(); ++i)
		{
			int tGlob, tLoc = vStepConv[i].numIter2a;
			ProcCom.allreduce(&tLoc, &tGlob, 1, PCL_DT_INT, PCL_RO_SUM);
			tGlob /= pcl::GetNumProcesses();

			if (!bPrintOnlyAverages) UG_LOG(std::setw(3) << tGlob << ",");

			sumAvgNumIter += tGlob;
			m_totalIterCntOfInnerSolvers += tGlob; // sum over all calls of inner solvers

			tLoc = vStepConv[i].numIter2a;
			ProcCom.allreduce(&tLoc, &tGlob, 1, PCL_DT_INT, PCL_RO_MAX);

			if (!bPrintOnlyAverages) UG_LOG(std::setw(3) << tGlob << ",");

			tLoc = vStepConv[i].numIter2a;
			ProcCom.allreduce(&tLoc, &tGlob, 1, PCL_DT_INT, PCL_RO_MIN);

			if (!bPrintOnlyAverages) UG_LOG(std::setw(3) << tGlob << "|");
		}
		if (bPrintOnlyAverages) {
			UG_LOG(" avg last defect: " << std::setprecision(2) << avgLastDef2a << " in ");
		}
		UG_LOG("" << std::setw(3) << sumAvgNumIter << " avg steps");
		if (type.compare("apply_F  ") == 0) {
//			UG_LOG(", " << std::fixed << (double)sumAvgNumIter/vStepConv.size()
//				        << std::scientific << " per FETI iteration");
			UG_LOG(" in " << std::setw(3) << vStepConv.size() << " call(s), "
				          << std::fixed << (double)sumAvgNumIter/vStepConv.size()
				          << std::scientific << " per FETI iteration");
		}
		UG_LOG("\n");
		
		if (bPrintOnlyAverages) UG_LOG("nps calls for '"<< type << "'          (7): ");

		if (!bPrintOnlyAverages) UG_LOG("Defect7   (avg)          :  ");
		for(size_t i = 0; i < vStepConv.size(); ++i)
		{
			double tGlob, tLoc = vStepConv[i].lastDef7;
			ProcCom.allreduce(&tLoc, &tGlob, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
			tGlob /= pcl::GetNumProcesses();
			avgLastDef7 = tGlob;

			if (!bPrintOnlyAverages) UG_LOG(std::setprecision(2) << tGlob << " |  ");
		}
		if (!bPrintOnlyAverages) {
			UG_LOG("\n");

			UG_LOG("NumIter7  (avg, max, min):");
		}
		sumAvgNumIter = 0;
		for(size_t i = 0; i < vStepConv.size(); ++i)
		{
			int tGlob, tLoc = vStepConv[i].numIter7;
			ProcCom.allreduce(&tLoc, &tGlob, 1, PCL_DT_INT, PCL_RO_SUM);
			tGlob /= pcl::GetNumProcesses();

			if (!bPrintOnlyAverages) UG_LOG(std::setw(3) << tGlob << ",");

			sumAvgNumIter += tGlob;
			m_totalIterCntOfInnerSolvers += tGlob; // sum over all calls of inner solvers

			tLoc = vStepConv[i].numIter7;
			ProcCom.allreduce(&tLoc, &tGlob, 1, PCL_DT_INT, PCL_RO_MAX);

			if (!bPrintOnlyAverages) UG_LOG(std::setw(3) << tGlob << ",");

			tLoc = vStepConv[i].numIter7;
			ProcCom.allreduce(&tLoc, &tGlob, 1, PCL_DT_INT, PCL_RO_MIN);

			if (!bPrintOnlyAverages) UG_LOG(std::setw(3) << tGlob << "|");
		}
		if (bPrintOnlyAverages) {
			UG_LOG(" avg last defect: " << std::setprecision(2) << avgLastDef7 << " in ");
		}
		UG_LOG("" << std::setw(3) << sumAvgNumIter << " avg steps");
		if (type.compare("apply_F  ") == 0) {
//			UG_LOG(", " << std::fixed << (double)sumAvgNumIter/vStepConv.size()
//				        << std::scientific << " per FETI iteration");
			UG_LOG(" in " << std::setw(3) << vStepConv.size() << " call(s), "
				          << std::fixed << (double)sumAvgNumIter/vStepConv.size()
				          << std::scientific << " per FETI iteration");
		}
		UG_LOG("\n");
} //
		
		if (!bPrintOnlyAverages) UG_LOG("\n");
	}
} /* end 'PrimalSubassembledMatrixInverse::print_statistic_of_inner_solver()' */


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	FETISolver implementation
template <typename TAlgebra>
FETISolver<TAlgebra>::
FETISolver() :
	m_spOperator(NULL),
	m_pMatrix(NULL),
	m_spDirichletSolver(NULL),
	m_spNeumannSolver(NULL)
{

}

template <typename TAlgebra>
bool FETISolver<TAlgebra>::
init(SmartPtr<MatrixOperator<matrix_type, vector_type> > A)
{
//	status
	UG_LOG("\n% Initializing FETI Solver: \n");

//	bool flag
	bool bSuccess = true;

//	remember A
	m_spOperator = A;

//	0. get matrix
	m_pMatrix = &m_spOperator->get_matrix();

//	check that DDInfo has been set
	if(m_pDDInfo == NULL)
	{
		UG_LOG("ERROR in FETISolver::init: DDInfo not set.\n");
		return false;
	}

	bool debugLayouts = (debug_writer()==NULL) ? false : true;

//	1. create FETI Layouts
	UG_LOG("\n%   - Create FETI layouts ... ");
	FETI_PROFILE_BEGIN(FETISolverInit_Create_Layouts);
	m_fetiLayouts.create_layouts(m_pMatrix->layouts(),
	                             m_pMatrix->num_rows(),
	                             *m_pDDInfo,
	                             debugLayouts);
	FETI_PROFILE_END();			// end 'FETI_PROFILE_BEGIN(FETISolverInit_Create_Layouts)' - Messpunkt ok
	UG_LOG("done.\n");

//  ----- 2. CONFIGURE LOCAL SCHUR COMPLEMENT  ----- //

//	2.1 set layouts in LocalSchurComplement
	m_LocalSchurComplement.set_feti_layouts(m_fetiLayouts);

//	2.2 init Dirichlet system and solver
//  check that dirichlet solver has been set
	if(m_spDirichletSolver.invalid())
	{
		UG_LOG("ERROR in FETISolver::init: No dirichlet solver set "
				" for inversion of A_{II} in Local Schur complement.\n");
		return false;
	}

//	set dirichlet solver for local Schur complement
	m_LocalSchurComplement.set_dirichlet_solver(m_spDirichletSolver);

//	set operator in local Schur complement
	m_LocalSchurComplement.set_matrix(m_spOperator);

//	2.3 init local Schur complement
	UG_LOG("\n%   - Init Local Schur Complement ... ");
	FETI_PROFILE_BEGIN(FETISolverInit_InitLocalSchurComplement);
	m_LocalSchurComplement.init();
	UG_LOG("done.\n");

//	2.4 check all procs
	if(!pcl::AllProcsTrue(bSuccess))
	{
		UG_LOG("ERROR in FETISolver::init: Some processes could not init"
				" local Schur complement.\n");
		return false;
	}
	FETI_PROFILE_END();			// end 'FETI_PROFILE_BEGIN(FETISolverInit_InitLocalSchurComplement)' - Messpunkt ok!

//  ----- 3. CONFIGURE SCHUR COMPLEMENT INVERSE  ----- //

//	3.1 set layouts in PrimalSubassembledMatrixInverse
	m_PrimalSubassembledMatrixInverse.set_feti_layouts(m_fetiLayouts);

//	3.2 init Neumann system and solver
//	check that neumann solver has been set
	if(m_spNeumannSolver.invalid())
	{
		UG_LOG("ERROR in FETISolver::init: No neumann solver set "
				" for inversion of A_{I,Delta}{I,Delta} in PrimalSubassembledMatrixInverse.\n");
		return false;
	}

//	set neumann solver in PrimalSubassembledMatrixInverse
	m_PrimalSubassembledMatrixInverse.set_neumann_solver(m_spNeumannSolver);

//  3.3 init coarse problem solver used in PrimalSubassembledMatrixInverse
//	check that coarse problem solver has been set
	if(m_spCoarseProblemSolver.invalid())
	{
		UG_LOG("ERROR in FETISolver::init: No coarse problem solver set "
				" for solving S_{Pi Pi} u_{Pi} = tilde{f}_{Pi}.\n");
		return false;
	}

//	set coarse problem solver in PrimalSubassembledMatrixInverse
	m_PrimalSubassembledMatrixInverse.set_coarse_problem_solver(m_spCoarseProblemSolver);

//	3.4 init PrimalSubassembledMatrixInverse (operator - given as parameter here - is also set thereby)
	FETI_PROFILE_BEGIN(FETISolverInit_InitPrimalSubassMatInv);
	UG_LOG("\n%   - Init 'PrimalSubassembledMatrixInverse' ...\n");
	m_PrimalSubassembledMatrixInverse.set_statistic_type("PSMI_init");
	if(m_PrimalSubassembledMatrixInverse.init(m_spOperator) != true)
	{
		UG_LOG("ERROR in FETISolver::init: Can not init Schur "
				"complement inverse.\n");
		bSuccess = false;
	}
	UG_LOG("\n%     'PrimalSubassembledMatrixInverse::init()' done.\n");

//	3.5 check all procs
	if(!pcl::AllProcsTrue(bSuccess))
	{
		UG_LOG("ERROR in FETISolver::init: Some processes could not init"
				" Schur complement inverse.\n");
		return false;
	}
	FETI_PROFILE_END();			// end 'FETI_PROFILE_BEGIN(FETISolverInit_InitPrimalSubassMatInv)' - Messpunkt ok!

//	status
	UG_LOG("\n% 'FETISolver::init()' done!\n");

//	we're done
	return true;
} /* end 'FETISolver::init()' */

template <typename TAlgebra>
bool FETISolver<TAlgebra>::
apply_return_defect(vector_type& u, vector_type& f)
{
//	status
	UG_LOG("\n% 'FETISolver::apply()':\n");
//	FETI_PROFILE_FUNC(); // should report same times as in section 'applyLinearSolver' (see 'operator_util.h')
//	FETI_PROFILE_BEGIN(FETISolverApplyReturnDefect); // profiling complete method
//	This function is used to solve the system Au=f. While the matrix A has
//	already been passed in additive storage in the function init(), here are
//	passed the (unknown, to be computed) solution u, that must be returned in
//	consistent storage, and the right-hand side, that is given in additive
//	storage. Note, that the right-hand side is overwritten, such that it can
//	be used for the computation of the defect (Once the defect is computed the
//	right-hand side is not needed anymore).

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

//	Lagrange multiplier
	vector_type lambda; lambda.create(u.size());
	m_fetiLayouts.vec_use_std_communication(lambda);

//	residuum
	vector_type r; r.create(u.size());
	m_fetiLayouts.vec_use_std_communication(r);

//	search direction
	vector_type p; p.create(u.size());
	m_fetiLayouts.vec_use_std_communication(p);
	//m_fetiLayouts.vec_use_???_communication(u); // unnec. since this setting is done before usage in 'm_PrimalSubassembledMatrixInverse.apply()' (see below)!

//	preconditioned residuum
	vector_type z; z.create(u.size());
	m_fetiLayouts.vec_use_std_communication(z);	// added 25022011ih

//	help vector to compute t = F*p
	vector_type t; t.create(u.size());
	m_fetiLayouts.vec_use_std_communication(t);

//	set start value of vector of Lagrange multipliers
	number lambdaStart = 0.0;
	bool isLambdaStartZero;
	if (lambdaStart == 0.0)	isLambdaStartZero = true; else 	isLambdaStartZero = false;
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

// 	Build start residuum:  r = r0 := d - F*lambda.
//	This is done in three steps.

//	(a) Compute d:= B S_{\Delta \Delta}^{-1} \tilde{f}_{\Delta}
//		and set r0 = d (prelimilarly, -F \lambda added later)
	FETI_PROFILE_BEGIN(FETISolverApply_Compute_D);
	m_PrimalSubassembledMatrixInverse.set_statistic_type("compute_d");
	if(!compute_d(r, f))
	{
		UG_LOG("ERROR in 'FETISolver::apply': "
			   "Cannot compute rhs 'd'. Aborting.\n");
		return false;
	}
	FETI_PROFILE_END();			// end 'FETI_PROFILE_BEGIN(FETISolverApply_Compute_D)' - Messpunkt ok, wenn 'ComputeDifferenceOnDelta()' als letzte Op. nicht stoert

	if (!isLambdaStartZero) {
// 	(b) Build t = F*lambda (t is additive afterwards)
		FETI_PROFILE_BEGIN(FETISolverApply_Apply_F);
		m_PrimalSubassembledMatrixInverse.set_statistic_type("apply_F  ");
		if(!apply_F(t, lambda))
		{
			UG_LOG("ERROR in 'FETISolver::apply': Unable "
				   "to build t = F*p. Aborting.\n");
			return false;
		}
		FETI_PROFILE_END();			// end 'FETI_PROFILE_BEGIN(FETISolverApply_Apply_F)' - Messpunkt ok, wenn 'ComputeDifferenceOnDelta()' als letzte Op. nicht stoert

// (c) Subtract values on \Delta, r0 = r0 - t
		m_fetiLayouts.vec_scale_append_on_dual(r, t, -1.0);
	}

//	prepare appearance of conv check
	prepare_conv_check();

//	compute and set start defect
	convergence_check()->start_defect(m_fetiLayouts.vec_norm_on_identified_dual(r));

// 	Precondition the start defect: apply z = M^-1 * r
	FETI_PROFILE_BEGIN(FETISolverApply_Apply_M_inverse);
	m_LocalSchurComplement.set_statistic_type("apply_M_inv_iter "); // maybe set type to "apply_M_inv_start" for debugging
	if (!apply_M_inverse(z, r)) // (calls 'm_LocalSchurComplement.apply()')
	{
		UG_LOG("ERROR in 'FETISolver::apply': "
			   "Cannot apply preconditioner. Aborting.\n");
		return false;
	}
	FETI_PROFILE_END();			// end 'FETI_PROFILE_BEGIN(FETISolverApply_Apply_M_inverse)' - Messpunkt ok, wenn 'ComputeDifferenceOnDelta()' + 'apply_scaling_matrix()' nicht stoeren

//	start values
	number rho, rho_new, beta, alpha, alpha_denominator;
	rho = rho_new = beta = alpha = alpha_denominator = 0.0;

// 	start search direction
	p = z;

// 	start rho
	rho = m_fetiLayouts.vec_prod_on_identified_dual(z, r);

//	very small number to check denominator
	number small = 1e-20;

// 	"lambda iteration" loop
	FETI_PROFILE_BEGIN(FETISolverApply_Lambda_iter_loop);
	while(!convergence_check()->iteration_ended())
	{
	//	increase iteration count
		m_iterCnt++;

	// 	Build t = F*p
		// p is consistent
		// t is consistent afterwards
		//m_fetiLayouts.vec_set_on_primal(p, 0.0); // seems to be not nec. (28022011ih)
		m_PrimalSubassembledMatrixInverse.set_statistic_type("apply_F  ");
		if(!apply_F(t, p))
		{
			UG_LOG("ERROR in 'FETISolver::apply': Unable "
						"to build t = F*p. Aborting.\n"); return false;
		}

	// 	Compute alpha
		alpha_denominator = m_fetiLayouts.vec_prod_on_identified_dual(t, p);

		if(fabs(alpha_denominator) < small)
		{
			UG_LOG("ERROR in 'FETISolver::apply': "
				   "VecProd(t, p) < small. Aborting.\n"); return false;
		}

		alpha = rho/alpha_denominator;

	// 	Update lambda := lambda + alpha*p
		m_fetiLayouts.vec_scale_add_on_dual(lambda, 1.0, lambda, alpha, p);

	// 	Update r := r - alpha*t
		m_fetiLayouts.vec_scale_add_on_dual(r, 1.0, r, -alpha, t);

	// 	Compute new defect
		convergence_check()->update_defect(m_fetiLayouts.vec_norm_on_identified_dual(r));
		if(convergence_check()->iteration_ended())
		{
			break;
			FETI_PROFILE_END();	// additional end 'FETI_PROFILE_BEGIN(FETISolverApply_Lambda_iter_loop)' - Messpunkt ok, da Konvergenz-Check ausgefuehrt
		}

	// 	Preconditioning: apply z = M^-1 * r
		m_LocalSchurComplement.set_statistic_type("apply_M_inv_iter ");
		if (!apply_M_inverse(z, r)) // (calls 'm_LocalSchurComplement.apply()')
		{
			UG_LOG("ERROR in 'FETISolver::apply': "
				   "Cannot apply preconditioner. Aborting.\n"); return false;
		}

	// 	new rho
		rho_new = m_fetiLayouts.vec_prod_on_identified_dual(z, r);

	// 	new beta
		beta = rho_new/rho;

	// 	new direction p:= beta*p + z
		m_fetiLayouts.vec_scale_add_on_dual(p, beta, p, 1.0, z);

	// 	update rho
		rho = rho_new;
	} /* end iteration loop */
	FETI_PROFILE_END();			// end 'FETI_PROFILE_BEGIN(FETISolverApply_Lambda_iter_loop)' - Messpunkt ok, da Konvergenz-Check ausgefuehrt

//	"back solve" (cf. A. Toselli, O. Widlund: "Domain Decomposition Methods -
//	Algorithms and Theory", chap. 6.4, p.164, l. 2) (03022011av)
	FETI_PROFILE_BEGIN(FETISolverApply_Backsolve);

//	reset t = 0.0
	t.set(0.0);

// compute t = B^T * lambda (Please note that 'ComputeDifferenceOnDeltaTransposed()'
// does not perform any communication!)
	ComputeDifferenceOnDeltaTransposed(t, lambda,
	                                   m_fetiLayouts.get_dual_master_indices(),
	                                   m_fetiLayouts.get_dual_slave_indices(),
	                                   m_fetiLayouts.get_dual_nbr_slave_indices());

	m_fetiLayouts.vec_use_intra_sd_communication(t);
	t.set_storage_type(PST_CONSISTENT);
	t.change_storage_type(PST_ADDITIVE);

//	compute f = f - B^T * lambda
	m_fetiLayouts.vec_scale_append_on_dual(f, t, -1.0);

	bool bSuccess = true;

//	Solve: A u = f
	FETI_PROFILE_BEGIN(FETISolverApply_ApplyPrimalSubassMatInv);
	m_PrimalSubassembledMatrixInverse.set_statistic_type("backsolve");
//	m_fetiLayouts.vec_use_intra_sd_communication(u); // added 25022011ih, but unnec. 01112011ih
//	m_fetiLayouts.vec_use_intra_sd_communication(f); // added 25022011ih, but unnec. 01112011ih
	if(!m_PrimalSubassembledMatrixInverse.apply(u, f))
	{
		UG_LOG("ERROR in FETISolver::apply: Cannot back solve.\n");
		bSuccess = false;
	}
	FETI_PROFILE_END();			// end 'FETI_PROFILE_BEGIN(FETISolverApply_ApplyPrimalSubassMatInv)' - Messpunkt ok, da 'AllProcsTrue()' am Ende von 'apply()' ausgefuehrt wurde

//	check all procs
	if(!pcl::AllProcsTrue(bSuccess))
	{
		UG_LOG("ERROR in FETISolver::apply: Some processes could not back solve.\n");
		return false;
	}
	FETI_PROFILE_END();			// end 'FETI_PROFILE_BEGIN(FETISolver_Backsolve)' - Messpunkt ok!

	return convergence_check()->post();
	//FETI_PROFILE_END();			// end 'FETI_PROFILE_BEGIN(FETISolverApplyReturnDefect)' - complete method

//	call this for output.
//	PROFILER_UPDATE();
//	PROFILER_OUTPUT("feti_profiling.rtf");
} /* end 'FETISolver::apply_return_defect()' */

template <typename TAlgebra>
bool FETISolver<TAlgebra>::
apply_F(vector_type& f, const vector_type& v)
{
	FETI_PROFILE_FUNC();
//	Help vector
	vector_type fTmp; fTmp.create(v.size());
	fTmp.set_storage_type(PST_CONSISTENT);

//	0. Reset values of f
	f.set(0.0);

//	1. Apply transposed jump operator: f = B_{\Delta}^T * v_{\Delta}:
//     v must be consistent (Please note that 'ComputeDifferenceOnDeltaTransposed()'
//     does not perform any communication!)
	ComputeDifferenceOnDeltaTransposed(f, v,
	                                   m_fetiLayouts.get_dual_master_indices(),
	                                   m_fetiLayouts.get_dual_slave_indices(),
	                                   m_fetiLayouts.get_dual_nbr_slave_indices());

//  2. Apply PrimalSubassembledMatrixInverse to f
	// f is consistent now, we make it additive
	m_fetiLayouts.vec_use_intra_sd_communication(f);   // nec. for 'change_storage_type()', but *not* for 'm_PrimalSubassembledMatrixInverse.apply()'!
	f.set_storage_type(PST_CONSISTENT);
	f.change_storage_type(PST_ADDITIVE);

	FETI_PROFILE_BEGIN(FETISolverApply_F_ApplyPrimalSubassMatInv);
	m_PrimalSubassembledMatrixInverse.apply(fTmp, f);
	FETI_PROFILE_END();			// end 'FETI_PROFILE_BEGIN(FETISolverApply_F_ApplyPrimalSubassMatInv)' - Messpunkt ok, da 'AllProcsTrue()' am Ende von 'apply()' ausgefuehrt wurde

//	3. Apply jump operator to get the final 'f'
	ComputeDifferenceOnDelta(f, fTmp, m_fetiLayouts.get_dual_master_layout(),
	                         	 	  m_fetiLayouts.get_dual_slave_layout(),
	                         	 	  m_fetiLayouts.get_dual_nbr_master_layout(),
	                         	 	  m_fetiLayouts.get_dual_nbr_slave_layout());

//	we're done
	return true;
} /* end 'FETISolver::apply_F()' */


template <typename TAlgebra>
bool FETISolver<TAlgebra>::
compute_d(vector_type& d, const vector_type& f)
{
	FETI_PROFILE_FUNC();
//	On entry, the vector f is filled with values for the dual unknowns. We make
//	no assumption on the values in the others (neglected unknowns).

//	create a help vector
	vector_type dTmp; dTmp.create(f.size());

//	0. Reset values of d, dTmp
	dTmp.set(0.0); d.set(0.0);

//  1. Apply PrimalSubassembledMatrixInverse to 'f'
//	1.1. let vectors use communication within feti subdomain - no 'vec_use_intra_sd_communication()' nec. before 'm_PrimalSubassembledMatrixInverse.apply()' (01112011ih)!

	dTmp.set_storage_type(PST_CONSISTENT);
	FETI_PROFILE_BEGIN(FETISolverCompute_d_ApplyPrimalSubassMatInv);
	if(!m_PrimalSubassembledMatrixInverse.apply(dTmp, f))
	{
		UG_LOG("In 'FETISolver::compute_d': Could not apply Schur"
				" complement inverse.\n");
		return false;
	}
	FETI_PROFILE_END();			// end 'FETI_PROFILE_BEGIN(FETISolverCompute_d_ApplyPrimalSubassMatInv)' - Messpunkt ok, da 'AllProcsTrue()' am Ende von 'apply()' ausgefuehrt wurde

//	2. Apply jump operator to get the final 'd'
	ComputeDifferenceOnDelta(d, dTmp, m_fetiLayouts.get_dual_master_layout(),
	                         	 	  m_fetiLayouts.get_dual_slave_layout(),
	                         	 	  m_fetiLayouts.get_dual_nbr_master_layout(),
	                         	 	  m_fetiLayouts.get_dual_nbr_slave_layout());

//	we're done
	return true;
} /* end 'FETISolver::compute_d()' */

template <typename TAlgebra>
bool FETISolver<TAlgebra>::
apply_M_inverse(vector_type& z, const vector_type& r)
{
	FETI_PROFILE_FUNC();
//	The incoming vector r is defined on the space V := range{B_{Delta}}. Thus,
//	on entry we assume the vector z to be consistent on the dual unknowns. The
//	primal and inner values are assumed to be undefined, since they are not needed
//	in the algorithm. Since z is consistent, in every dualMaster<->dualSlave is
//	stored the same value and in every dualNbrMaster<->dualNbrSlave is stored
//	the same value.

//	Help vector
	vector_type zTmp; zTmp.create(r.size()); zTmp = z;

//	0. Reset values of z, zTmp
	z.set(0.0); zTmp.set(0.0);

//	1. Apply scaling: z := D_{\Delta}^{(i)} * r
	apply_scaling_matrix(z, r); // maybe restrict to layout (but it performs no communication until now)

//  2. 	Apply transposed jump operator: zTmp := B_{\Delta}^T * z
//		Afterwards, zTmp is consistent on Delta (Please note that 'ComputeDifferenceOnDeltaTransposed()'
//      does not perform any communication!)
	ComputeDifferenceOnDeltaTransposed(zTmp, z,
	                                   m_fetiLayouts.get_dual_master_indices(),
	                                   m_fetiLayouts.get_dual_slave_indices(),
	                                   m_fetiLayouts.get_dual_nbr_slave_indices());

//	3. Apply local Schur complement: z := S_{\Delta}^{(i)} * zTmp
//	3.1. let vectors use communication within feti subdomain
	m_fetiLayouts.vec_use_intra_sd_communication(z);
	m_fetiLayouts.vec_use_intra_sd_communication(zTmp);
//	3.2. set correct parallel storage type
	z.set_storage_type(PST_ADDITIVE);
	zTmp.set_storage_type(PST_CONSISTENT);
//	3.3. solve
	FETI_PROFILE_BEGIN(FETISolverApply_M_inv_ApplyLocalSchurComplement);
	m_LocalSchurComplement.apply(z, zTmp);
	FETI_PROFILE_END();			// end 'FETI_PROFILE_BEGIN(FETISolverApply_M_inv_ApplyLocalSchurComplement)' - Messpunkt ok, da 'AllProcsTrue()' am Ende von 'apply()' aufgerufen wird.
								// ('m_pOperator->apply' allein synchronisiert *nicht*!)
//  4. Apply jump operator:  zTmp :=  B_{\Delta} * z
	ComputeDifferenceOnDelta(zTmp, z, m_fetiLayouts.get_dual_master_layout(),
	                         	 	  m_fetiLayouts.get_dual_slave_layout(),
	                         	 	  m_fetiLayouts.get_dual_nbr_master_layout(),
	                         	 	  m_fetiLayouts.get_dual_nbr_slave_layout());

//	5. Apply scaling: z := D_{\Delta}^{(i)} * zTmp to get the final 'z'
	apply_scaling_matrix(z, zTmp); // maybe restrict to layout (but it performs no communication until now)

//	we're done
	return true;
} /* end 'FETISolver::apply_M_inverse()' */

template <typename TAlgebra>
void FETISolver<TAlgebra>::
test_layouts(bool bPrint)
{
	m_fetiLayouts.test_layouts(bPrint);
} /* end 'FETISolver::test_layouts()' */

////////////////////////////////////////////////////////////////////////
//	template instantiations for all current algebra types.




#ifdef UG_CPU_1
template class LocalSchurComplement<CPUAlgebra>;
template class PrimalSubassembledMatrixInverse<CPUAlgebra>;
template class FETISolver<CPUAlgebra>;
#endif
#ifdef UG_CPU_2
template class LocalSchurComplement<CPUBlockAlgebra<2> >;
template class PrimalSubassembledMatrixInverse<CPUBlockAlgebra<2> >;
template class FETISolver<CPUBlockAlgebra<2> >;	
#endif
#ifdef UG_CPU_3
template class LocalSchurComplement<CPUBlockAlgebra<3> >;
template class PrimalSubassembledMatrixInverse<CPUBlockAlgebra<3> >;
template class FETISolver<CPUBlockAlgebra<3> >;
#endif
#ifdef UG_CPU_4
template class LocalSchurComplement<CPUBlockAlgebra<4> >;
template class PrimalSubassembledMatrixInverse<CPUBlockAlgebra<4> >;
template class FETISolver<CPUBlockAlgebra<4> >;
#endif
#ifdef UG_CPU_VAR
template class LocalSchurComplement<CPUVariableBlockAlgebra>;
template class PrimalSubassembledMatrixInverse<CPUVariableBlockAlgebra >;
template class FETISolver<CPUVariableBlockAlgebra >;
#endif
};  // end of namespace

#endif
