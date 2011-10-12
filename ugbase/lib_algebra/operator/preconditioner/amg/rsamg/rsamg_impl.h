/**
 * \file rsamg_impl.h
 *
 * \author Martin Rupp
 *
 * \date 03.12.2009
 *
 * implementation file for amg
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#ifndef __H__UG__LIB_ALGEBRA__RSAMG_SOLVER__RSAMG_IMPL_H__
#define __H__UG__LIB_ALGEBRA__RSAMG_SOLVER__RSAMG_IMPL_H__

//#include "sparsematrix_util.h"

#include "rsamg.h"
#include "../graph.h"

#include "rsamg_prolongation.h"
#include "rsamg_debug.h"
#include "../amg_debug_helper.h"
#include "rsamg_nodeinfo.h"

#include "rsamg_coarsening.h"

#include "../stopwatch.h"
#include "../amg_profiling.h"

#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#include "lib_algebra/parallelization/parallelization_util.h"
#include "rsamg_parallel_coarsening.h"
#endif

std::string GetProcFilename(std::string path, std::string name, std::string extension);


namespace ug{

//#define GRAPH_WITH_LOCAL_INVERSE


#define AMG_WRITE_MATRICES_PATH "/Users/mrupp/matrices/AMG_"
#define AMG_WRITE_MATRICES_MAX (200*200)



#if 0
#define AMG_WRITE_GRAPH


//#define AMG_WRITE_GRAPH


//#define AMG_PRINT_INDIRECT

#define AMG_PRINT_GRAPH

#define AMG_WRITE_COARSENING

#define AMG_PRINT_COARSENING
#define AMG_PRINT_P
#define AMG_PRINT_R
#define AMG_PRINT_AH

#define AMG_PRINT_COARSEN_RATINGS
#define AMG_PRINT_COARSEN
#endif

inline double amg_diag_value(const double &d) { return d; }
inline double amg_offdiag_value(const double &d) { return d; }

template<typename T> inline double amg_diag_value(const T &d) { return BlockNorm(d); }
template<typename T> inline double amg_offdiag_value(const T &d) { return -BlockNorm(d); }


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateStrongConnectionGraph:
//------------------------------
/**
 * Create graph of strong connections from a matrix A
 * i is strong coupled to j, if \f$ -A_{ij} \ge m_dTheta {\max\limits}_{\phi_1(A_{ik}) < 0} { \phi_2(A_{ik}) }\f$
 * at the moment, we use the function \f$\phi_1 = \phi_2 =\f$ amg_offdiag_value
 * \param	A			matrix A for which to calculate strong connectivity graph
 * \param 	graph		the calculated strong connectivity graph of A
 *						graph is afterwards made up of connections from a node i to j if
 *						j has a strong connection to i
 * \param	m_dTheta	the smaller, the more connections are considered strong e.g. 0.25 is a reasonable default value
 * \note	you can also pre-calculate some connectivity-measurement matrix C, and call this function with C.
 */
template<typename matrix_type>
void CreateStrongConnectionGraph(const matrix_type &A, cgraph &graph, double m_dTheta)
{
	AMG_PROFILE_FUNC();
	graph.resize(A.num_rows());

	for(size_t i=0; i< A.num_rows(); i++)
	{
		if(A.is_isolated(i))
			continue;

		double dmax = 0;

		for(typename matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
		{
			if(conn.index() == i) continue; // skip diag
			if(conn.value() != 0.0 && amg_offdiag_value(conn.value()) < dmax)
				dmax = amg_offdiag_value(conn.value());
		}

		for(typename matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
		{
			if(conn.index() == i) continue; // skip diag
			if( amg_offdiag_value(conn.value()) < m_dTheta * dmax)
				graph.set_connection(i, conn.index());
		}

		UG_ASSERT(graph.num_connections(i) > 0, "");
	}

#ifdef AMG_PRINT_GRAPH
	graph.print();
#endif
}


/*template<typename TAlgebra>
void rsamg<TAlgebra>::FullSubdomainBlocking(matrix_type &AH, IndexLayout &nextMasterLevelLayout, IndexLayout &nextSlaveMasterLayout)
{

}*/

template<typename TAlgebra>
void RSAMG<TAlgebra>::create_new_indices(stdvector<int> &newIndex, const AMGNodes &nodes, size_t level)
{
	AMG_PROFILE_FUNC();
	std::vector<bool> &vIsFine = levels[level]->is_fine;
	vIsFine.resize(nodes.size());
	size_t iNrOfCoarse=0;
	for(size_t i=0; i < nodes.size(); i++)
		if(nodes[i].is_coarse())
		{
			newIndex[i] = iNrOfCoarse;
			vIsFine[i] = false;
			iNrOfCoarse++;
		}
		else
		{
			vIsFine[i] = true;
			newIndex[i] = -1;
		}
	UG_ASSERT(iNrOfCoarse == nodes.get_nr_of_coarse(), "");
}

template<typename TAlgebra>
void RSAMG<TAlgebra>::create_parentIndex(const stdvector<int> &newIndex, const AMGNodes &nodes, size_t level)
{
	AMG_PROFILE_FUNC();
	m_parentIndex.resize(level+2);
	m_parentIndex[level+1].resize(nodes.get_nr_of_coarse());
	for(size_t i=0; i<nodes.size(); i++)
		if(nodes[i].is_coarse())
			m_parentIndex[level+1][ newIndex[i] ] = i;
}



template<typename TAlgebra>
template<typename TMatrix>

void RSAMG<TAlgebra>::write_debug_matrix(TMatrix &mat, size_t fromlevel, size_t tolevel, const char *name)
{
	if(m_writeMatrices)
	{
		AMG_PROFILE_FUNC();
		std::string filename = GetProcFilename(m_writeMatrixPath, ToString(name) + ToString(fromlevel),".mat");
		AMGWriteToFile(mat, fromlevel, tolevel, filename.c_str(), m_amghelper);
		std::fstream f2(filename.c_str(), std::ios::out | std::ios::app);
		f2 << "c " << GetProcFilename("", std::string("AMG_fine_L") + ToString(fromlevel), ".marks") << "\n";
		//f2 << "c " << GetProcFilename("", std::string("AMG_aggfine_L") + ToString(level), ".marks") << "\n";
		f2 << "c " << GetProcFilename("", std::string("AMG_coarse_L") + ToString(fromlevel), ".marks") << "\n";
		f2 << "c " << GetProcFilename("", std::string("AMG_other_L") + ToString(fromlevel), ".marks") << "\n";
		f2 << "c " << GetProcFilename("", std::string("AMG_dirichlet_L") + ToString(fromlevel), ".marks") << "\n";
		//f2 << "v " << GetProcFilename("", std::string("AMG_d_L") + ToString(fromlevel), ".marks") << "\n";
	}
}



template<typename TAlgebra>
void RSAMG<TAlgebra>::write_debug_matrix_markers
	(size_t level, const AMGNodes &nodes)
{
	AMG_PROFILE_FUNC();
	std::fstream ffine(GetProcFilename(m_writeMatrixPath, std::string("AMG_fine_L") + ToString(level), ".marks").c_str(), std::ios::out);
	ffine << "1 0 0 1 0\n";
	std::fstream ffine2(GetProcFilename(m_writeMatrixPath, std::string("AMG_aggfine_L") + ToString(level), ".marks").c_str(), std::ios::out);
	ffine2 << "1 0.2 1 1 0\n";
	std::fstream fcoarse(GetProcFilename(m_writeMatrixPath, std::string("AMG_coarse_L") + ToString(level), ".marks").c_str(), std::ios::out);
	fcoarse << "0 0 1 1 2\n";
	std::fstream fother(GetProcFilename(m_writeMatrixPath, std::string("AMG_other_L") + ToString(level), ".marks").c_str(), std::ios::out);
	fother << "1 1 0 1 0\n";
	std::fstream fdirichlet(GetProcFilename(m_writeMatrixPath, std::string("AMG_dirichlet_L") + ToString(level), ".marks").c_str(), std::ios::out);
	fdirichlet << "0 1 1 1 0\n";
	/*for(size_t i=0; i < rating.size(); i++)
	{
		//int o = m_amghelper.GetOriginalIndex(level, i);
		int o=i;
		if(rating[i].is_aggressive_fine()) ffine2 << o << "\n";
		else if(rating[i].is_fine()) ffine << o << "\n";
		else if(rating[i].is_coarse()) fcoarse << o << "\n";
		else if(rating[i].is_dirichlet()) fdirichlet << o << "\n";
		else fother << o << "\n";
	}*/
	for(size_t i=0; i<nodes.size(); i++)
	{
		int o = i; //m_amghelper.GetOriginalIndex(level, i);
		if(nodes[i].is_fine_direct()) ffine << o << "\n";
		else if(nodes[i].is_coarse()) fcoarse << o << "\n";
		else fother << o << "\n";
	}
}

template<typename TAlgebra>
void RSAMG<TAlgebra>::debug_matrix_write(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
		prolongation_matrix_type &P, size_t level, const AMGNodes &nodes)
{
	if(m_amghelper.positions[level].size() > 0)
	{
		m_amghelper.positions.resize(level+2);
		m_amghelper.positions[level+1].resize(AH.num_rows());
		for(size_t i=0; i < AH.num_rows(); i++)
			m_amghelper.positions[level+1][i] = m_amghelper.positions[level][m_parentIndex[level+1][i]];
	}

	write_debug_matrix_markers(level, nodes);

	AMG_PROFILE_FUNC();
	UG_LOG("write matrices");

	write_debug_matrix(A, level, level, "AMG_A_L");				UG_LOG(".");
	write_debug_matrix(P, level, level, "AMG_P_L");	UG_LOG(".");
	write_debug_matrix(R, level, level+1, "AMG_R_L");			UG_LOG(".");

	AMGWriteToFile(AH, level+1, level+1, GetProcFilename(m_writeMatrixPath,
			ToString("AMG_A_L") + ToString(level+1),".mat").c_str(), m_amghelper);

	UG_LOG(". done\n");

	AMGWriteToFile(AH, level+1, level+1, GetProcFilename(m_writeMatrixPath,
					ToString("AMG_A_L") + ToString(level+1),".mat").c_str(), m_amghelper);

	UG_DLOG(LIB_ALG_AMG, 1, ". done.\n");
}


#ifdef UG_PARALLEL
template<typename matrix_type>
bool _MakeFullRowsMatrix(const ParallelMatrix<matrix_type> &_mat, ParallelMatrix<matrix_type> &newMat,
		IndexLayout &OL1MasterLayout, IndexLayout &OL1SlaveLayout)
{
	PROFILE_FUNC();
	// pcl does not use const much
	//UG_ASSERT(overlap_depth > 0, "overlap_depth has to be > 0");
	ParallelMatrix<matrix_type> &mat = const_cast<ParallelMatrix<matrix_type> &> (_mat);

	IndexLayout totalMasterLayout, totalSlaveLayout;
	stdvector<IndexLayout> vMasterLayouts, vSlaveLayouts;
	ParallelNodes PN(mat.get_communicator(), mat.get_master_layout(), mat.get_slave_layout(), mat.num_rows());
	GenerateOverlapClass<ParallelMatrix<matrix_type> > c(mat, newMat,
			totalMasterLayout, totalSlaveLayout, vMasterLayouts, vSlaveLayouts, PN);
	c.m_overlapDepthMaster = 1;
	c.m_overlapDepthSlave = 0;
	c.m_masterDirichletLast = false;
	c.m_slaveDirichletLast = false;
	bool b = c.calculate();

	AddLayout(OL1MasterLayout, mat.get_master_layout());
	AddLayout(OL1SlaveLayout, mat.get_slave_layout());
	AddLayout(OL1MasterLayout, vMasterLayouts[0]);
	AddLayout(OL1SlaveLayout, vSlaveLayouts[0]);
	//overlapSize = c.m_overlapSize;
	return b;
}
#endif


#ifdef UG_PARALLEL
void CreateAllToAllFromMasterSlave(pcl::ParallelCommunicator<IndexLayout> &communicator,
		IndexLayout &OLCoarseningSendLayout, IndexLayout &OLCoarseningReceiveLayout,
		IndexLayout &OL1MasterLayout, IndexLayout &OL1SlaveLayout);

#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// createAMGLevel:
//-------------------------
/**
 * create AMG matrix R, P, and AH = R A P
 * \param AH
 * \param R
 * \param A
 * \param P
 * \param level
 */
template<typename TAlgebra>
void RSAMG<TAlgebra>::create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
		prolongation_matrix_type &P, size_t level)
{
	AMG_PROFILE_FUNC();
#ifdef UG_PARALLEL
	//UG_ASSERT(pcl::GetNumProcesses()==1, "not implemented for procs > 1");


	matrix_type AOL1;
	IndexLayout OL1MasterLayout, OL1SlaveLayout;
	_MakeFullRowsMatrix(A, AOL1, OL1MasterLayout, OL1SlaveLayout);

	stdvector<bool> bMaster(AOL1.num_rows(), false);
	for(IndexLayout::iterator iter = A.get_master_layout().begin(); iter != A.get_master_layout().end(); ++iter)
	{
		IndexLayout::Interface &interface = A.get_master_layout().interface(iter);
		for(IndexLayout::Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
			bMaster[interface.get_element(iter2)] = true;
	}

	stdvector<bool> bSlave(AOL1.num_rows(), false);
	for(IndexLayout::iterator iter = A.get_slave_layout().begin(); iter != A.get_slave_layout().end(); ++iter)
	{
		IndexLayout::Interface &interface = A.get_master_layout().interface(iter);
		for(IndexLayout::Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
			bMaster[interface.get_element(iter2)] = true;
	}

#else
	const matrix_type &AOL1 = A;
#endif



	size_t NOL1 = AOL1.num_rows();
	size_t N = A.num_rows();
	AMGNodes nodes(NOL1);
	UG_LOG("NOL1 = " << NOL1 << "\n");

#ifdef UG_PARALLEL
	nodes.bMaster = bMaster;
	nodes.bSlave = bSlave;
#endif

	bool bTiming=true;
	UG_DLOG(LIB_ALG_AMG, 1, "Creating level " << level << ". (" << N << " nodes, " << NOL1-N << " overlapping )" << std::endl << std::fixed);
	stopwatch SW;

	// todo: check for isolated condition

	nodeinfo_pq_type PQ;

	stdvector<int> posInConnections; posInConnections.resize(NOL1, -1);

	UG_DLOG(LIB_ALG_AMG, 1, "A.totalNrOfConnections = " << A.total_num_connections() << std::endl);

	cgraph graphS(NOL1);
	cgraph graphST;

	// build graph
	/////////////////////////////////////////

	UG_DLOG(LIB_ALG_AMG, 1, "building graph... ");
	if(bTiming) SW.start();

	CreateStrongConnectionGraph(AOL1, graphS);

	// we need the transpose, since when we set a node coarse, we want
	// all nodes to be fine which can be interpolated by this coarse node
	// graph is afterwards made up of connections from a node i to j if
	// j has a strong connection to i
	graphST.set_as_transpose_of(graphS);

#ifdef AMG_WRITE_GRAPH
	WriteAMGGraphToFile(graphS, (std::string(AMG_WRITE_MATRICES_PATH) + "G" + ToString(level) + ".mat").c_str(), m_amghelper, level);
	WriteAMGGraphToFile(graphST, (std::string(AMG_WRITE_MATRICES_PATH) + "GT" + ToString(level) + ".mat").c_str(), m_amghelper, level);
#endif

	CreateMeasureOfImportancePQ(graphS, graphST, PQ, nodes);

	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

	IF_DEBUG(LIB_ALG_AMG, 4)
		nodes.print_ratings(m_amghelper, level);

	// Coarsen
	/////////////////////////////////////////
	UG_DLOG(LIB_ALG_AMG, 1, std::endl << "coarsening... ");

	SW.start();

#ifdef UG_PARALLEL
	// subdomain blocking
	enum ParallelCoarseningType
	{
		PCT_FULL_SUBDOMAIN_BLOCKING, PCT_COLORING_COARSEN
	};
	ParallelCoarseningType m_parallelCoarseningType = PCT_FULL_SUBDOMAIN_BLOCKING;

	if(m_parallelCoarseningType == PCT_FULL_SUBDOMAIN_BLOCKING)
	{
		FullSubdomainBlocking(graphST, PQ, nodes);
	}
	else if(m_parallelCoarseningType == PCT_COLORING_COARSEN)
	{
		IndexLayout OLCoarseningSendLayout, OLCoarseningReceiveLayout;

		CreateAllToAllFromMasterSlave(AOL1.get_communicator(), OLCoarseningSendLayout, OLCoarseningReceiveLayout, OL1MasterLayout, OL1SlaveLayout);

		ColoringCoarsen(AOL1.get_communicator(), OLCoarseningSendLayout, OLCoarseningReceiveLayout,
			graphST, PQ, nodes);
	}
#else
	Coarsen(graphST, PQ, nodes);
#endif

	PreventFFConnections(graphS, graphST, nodes);

	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

	// agressive Coarsening
	/////////////////////////////////////////

	if(m_bAggressiveCoarsening && level == 0)
	{
		// build graph 2
		//------------------

		UG_DLOG(LIB_ALG_AMG, 1, std::endl << "building graph2... ");
		if(bTiming) SW.start();

		UG_DLOG(LIB_ALG_AMG, 1, "unassigned = " << nodes.get_unassigned() << "\n");

		//unassigned = 0; ??
		cgraph graphAC(N);
		CreateAggressiveCoarseningGraph(graphST, graphAC, nodes, m_iAggressiveCoarseningNrOfPaths, &posInConnections[0]);
		CreateMeasureOfImportanceAggressiveCoarseningPQ(graphAC, PQ, nodes);
		if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");
		// coarsen 2
		//------------------

		if(nodes.get_unassigned() == 0)
		{
			UG_DLOG(LIB_ALG_AMG, 1, std::endl << "skipping coarsening2: no unassigned nodes.");
		}
		else
		{
			UG_DLOG(LIB_ALG_AMG, 1, std::endl << "coarsening2... ");

			if(bTiming) SW.start();
			Coarsen(graphAC, PQ, nodes);
			//PreventFFConnections(graphS, graphST, nodes);
			if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

		}
	}

	// todo: perhaps PreventFFConnections should be done HERE, after AC.

#ifdef AMG_PRINT_COARSENING
	printCoarsening(level);
#endif

	// construct prolongation P = I_{2h->h}
	/////////////////////////////////////////

	UG_DLOG(LIB_ALG_AMG, 1, std::endl << "prolongation... ");


	if(bTiming) SW.start();
	CreateRugeStuebenProlongation(P, A, nodes, //newIndex,
			m_dTheta, m_dEpsilonTr);
	UG_ASSERT(nodes.get_unassigned() == 0 || (m_bAggressiveCoarsening && level == 0),
			"no aggressive coarsening, but indirect interpolation of " << nodes.get_unassigned() << " nodes ?");
	if(nodes.get_unassigned() > 0)
		CreateIndirectProlongation(P, A, nodes, &posInConnections[0], m_dTheta);

	//----------

	/// create layouts


	stdvector<int> newIndex(nodes.size());
	create_new_indices(newIndex, nodes, level);

	// create parentIndex (for debug only)
	if(m_writeMatrices)
		create_parentIndex(newIndex, nodes, level);


	P.defragment();
#ifdef UG_PARALLEL
	P.set_storage_type(PST_CONSISTENT);
#endif

	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

#ifdef AMG_PRINT_P
	UG_DLOG(LIB_ALG_AMG, 1, std::endl << "Prolongation level " << level << std::endl);
	P.print();
#endif


	// construct restriction R = I_{h->2h}
	/////////////////////////////////////////

	UG_DLOG(LIB_ALG_AMG, 1, std::endl << "restriction... ");
	if(bTiming) SW.start();

	PROFILE_BEGIN(AMGSetAsTransposeOf)
	// construct restriction R = I_{h -> 2h}
	R.set_as_transpose_of(P);
	PROFILE_END();
	// R is already defragmented

#ifdef UG_PARALLEL
	R.set_storage_type(PST_CONSISTENT);
#endif

	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

#ifdef AMG_PRINT_R
	UG_DLOG(LIB_ALG_AMG, 1, std::endl << "Restriction level " << level << std::endl);
	R.print();
#endif

	// create Galerkin product
	/////////////////////////////////////////

	UG_DLOG(LIB_ALG_AMG, 1, "\ngalerkin product... ");
	if(bTiming) SW.start();

	// AH = R A P
	PROFILE_BEGIN(AMGCreateAsMultiplyOf)
		CreateAsMultiplyOf(AH, R, A, P);
	PROFILE_END();

	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

	// defragment
	if(bTiming) SW.start();
	AH.defragment();
	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, std::endl << "Defragment ... took " << SW.ms() << " ms");

#ifdef UG_PARALLEL
	AH.set_storage_type(PST_ADDITIVE);
	AH.set_layouts(A.get_master_layout(), A.get_slave_layout());
	AH.set_communicator((const_cast<matrix_type&>(A)).get_communicator());
	AH.set_process_communicator((const_cast<matrix_type&>(A)).get_process_communicator());
#endif

#ifdef AMG_PRINT_AH
	UG_DLOG(LIB_ALG_AMG, 1, "AH level " << level << std::endl);
	AH.print();
#endif

	UG_DLOG(LIB_ALG_AMG, 1, "\n");

	if(m_writeMatrices && A.num_rows() < AMG_WRITE_MATRICES_MAX)
		debug_matrix_write(AH, R, AOL1, P, level, nodes);
}


//!
//! amg constructor
template<typename TAlgebra>
RSAMG<TAlgebra>::RSAMG() :
	AMGBase<TAlgebra>(),
	m_dEpsilonTr(0.3),
	m_dTheta(0.3),
	//m_dSigma(0.3),
	m_bAggressiveCoarsening(0),
	m_iAggressiveCoarseningNrOfPaths(2)
{
}



template<typename TAlgebra>
void RSAMG<TAlgebra>::tostring() const
{
	AMGBase<TAlgebra>::tostring();

	UG_LOG("Ruge/Stueben AMG Preconditioner:\n");
	UG_LOG(" theta (epsilon_strong, strong connectivity) = " << m_dTheta << std::endl);
	//UG_LOG(" sigma = " << m_dSigma << std::endl);
	UG_LOG(" epsilon_tr (truncation of interpolation) = " << m_dEpsilonTr << std::endl);

	if(m_bAggressiveCoarsening)	{UG_LOG(" Aggressive Coarsening is on, A" << m_iAggressiveCoarseningNrOfPaths << "-mode." << std::endl);}
	else						{UG_LOG(" no Aggressive Coarsening" << std::endl);}
}

} // namespace ug

#endif //  __H__UG__LIB_ALGEBRA__RSAMG_SOLVER__RSAMG_IMPL_H__
