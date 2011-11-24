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


inline double amg_diag_value(const double &d) { return d; }
inline double amg_offdiag_value(const double &d) { return d; }

template<typename T> inline double amg_diag_value(const T &d) { return BlockNorm(d); }
template<typename T> inline double amg_offdiag_value(const T &d) { return -BlockNorm(d); }

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
#include "lib_algebra/parallelization/parallel_matrix_overlap_impl.h"
#include "rsamg_parallel_coarsening.h"
#endif


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


#ifdef UG_PARALLEL
template<typename matrix_type>
bool _MakeFullRowsMatrix(const ParallelMatrix<matrix_type> &_mat, ParallelMatrix<matrix_type> &newMat,
		ParallelNodes &PN, IParallelCoarsening *pParallelCoarsening,
		stdvector<IndexLayout> &vMasterLayouts, stdvector<IndexLayout> &vSlaveLayouts)
{
	PROFILE_FUNC();
	// pcl does not use const much
	//UG_ASSERT(overlap_depth > 0, "overlap_depth has to be > 0");
	ParallelMatrix<matrix_type> &mat = const_cast<ParallelMatrix<matrix_type> &> (_mat);

	IndexLayout totalMasterLayout, totalSlaveLayout;

	GenerateOverlapClass<ParallelMatrix<matrix_type> > c(mat, newMat,
			totalMasterLayout, totalSlaveLayout, vMasterLayouts, vSlaveLayouts, PN);
	c.m_overlapDepthMaster = pParallelCoarsening->overlap_depth_master();
	c.m_overlapDepthSlave = pParallelCoarsening->overlap_depth_slave();
	c.m_masterDirichletLast = false;
	c.m_slaveDirichletLast = false;
	bool b = c.calculate();

	//PRINTLAYOUT(mat.get_communicator(), mat.get_master_layout(), mat.get_slave_layout());
	//PRINTLAYOUT(mat.get_communicator(), vMasterLayouts[0], vSlaveLayouts[0]);

	//overlapSize = c.m_overlapSize;
	return b;
}
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
		prolongation_matrix_type &PnewIndices, size_t level)
{
	AMG_PROFILE_FUNC();

#ifdef UG_PARALLEL
	//UG_ASSERT(pcl::GetNumProcesses()==1, "not implemented for procs > 1");


	matrix_type AOL1;
	matrix_type &A_ = const_cast<matrix_type&> (A);
	ParallelNodes PN(A_.get_communicator(), A_.get_master_layout(), A_.get_slave_layout(), A_.num_rows());
	UG_ASSERT(m_pParallelCoarsening!=NULL, "please set the Parallel Coarsening type");

	stdvector<IndexLayout> vMasterLayouts, vSlaveLayouts;
	_MakeFullRowsMatrix(A, AOL1, PN, m_pParallelCoarsening, vMasterLayouts, vSlaveLayouts);


#else
	const matrix_type &AOL1 = A;
#endif
	prolongation_matrix_type PoldIndices;



	size_t NOL1 = AOL1.num_rows();
	size_t N = A.num_rows();
#ifdef UG_PARALLEL
	AMGNodes nodes(NOL1, PN);
#else
	AMGNodes nodes(NOL1);
#endif
	UG_LOG("NOL1 = " << NOL1 << "\n");

	bool bTiming=true;
	//UG_DLOG(LIB_ALG_AMG, 0, "Creating level " << level << ". (" << N << " nodes, " << NOL1-N << " overlapping )" << std::endl << std::fixed);
	UG_LOG("Creating level " << level << ". (" << N << " nodes, " << NOL1-N << " overlapping )" << std::endl << std::fixed);
	stopwatch SW;

	// todo: check for isolated condition

	nodeinfo_pq_type PQ;

	//A.print();

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

	bool bUnsymmetric = true;
#ifdef UG_PARALLEL
	m_pParallelCoarsening->coarsen(PN, vMasterLayouts, vSlaveLayouts,
			graphS, graphST, PQ, nodes, bUnsymmetric);
#else
	Coarsen(graphS, graphST, PQ, nodes, bUnsymmetric, false);
#endif

	//PreventFFConnections(graphS, graphST, nodes);
	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

	// agressive Coarsening
	/////////////////////////////////////////

	if(m_bAggressiveCoarsening && level == 0)
	{
#ifdef UG_PARALLEL
		UG_ASSERT(pcl::GetNumProcesses() == 1, "");
#endif
		// build graph 2
		//------------------

		UG_DLOG(LIB_ALG_AMG, 1, std::endl << "building graph2... ");
		if(bTiming) SW.start();

		UG_DLOG(LIB_ALG_AMG, 1, "unassigned = " << nodes.get_unassigned() << "\n");

		// todo: unsymmetric case.
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
			Coarsen(graphAC, graphAC, PQ, nodes, bUnsymmetric, true);
			//PreventFFConnections(graphS, graphST, nodes);
			if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

		}
	}

	// todo: perhaps PreventFFConnections should be done HERE, after AC.

	//nodes.print();
#ifdef AMG_PRINT_COARSENING
	printCoarsening(level);
#endif

	// construct prolongation P = I_{2h->h}
	/////////////////////////////////////////

	UG_DLOG(LIB_ALG_AMG, 1, std::endl << "prolongation... ");


	if(bTiming) SW.start();
	CreateRugeStuebenProlongation(PoldIndices, AOL1, nodes, //newIndex,
			m_dTheta, m_dEpsilonTr);
	if(!(nodes.get_unassigned() == 0 || (m_bAggressiveCoarsening && level == 0)))
	{
		nodes.print_ratings(m_amghelper, level);
		UG_ASSERT(nodes.get_unassigned() == 0 || (m_bAggressiveCoarsening && level == 0),"no aggressive coarsening, but indirect interpolation of " << nodes.get_unassigned() << " nodes ?");
	}

	if(nodes.get_unassigned_indirect_fine() > 0)
		CreateIndirectProlongation(PoldIndices, AOL1, nodes, &posInConnections[0], m_dTheta);

	//----------

	/// create layouts

#ifdef UG_PARALLEL
	this->parallel_process_prolongation(PoldIndices, PnewIndices, m_dEpsilonTr, level, nodes,
			PN, true, AH.get_master_layout(), AH.get_slave_layout());
#else
	this->serial_process_prolongation(PoldIndices, PnewIndices, m_dEpsilonTr, level, nodes);
#endif

	// construct restriction R = I_{h->2h}
	/////////////////////////////////////////

	UG_DLOG(LIB_ALG_AMG, 1, std::endl << "restriction... ");
	if(bTiming) SW.start();

	PROFILE_BEGIN(AMGSetAsTransposeOf)
	// construct restriction R = I_{h -> 2h}
	R.set_as_transpose_of(PnewIndices);
	PROFILE_END();
	// R is already defragmented

#ifdef UG_PARALLEL
	R.set_storage_type(PST_CONSISTENT);
#endif

	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");


	// create Galerkin product
	/////////////////////////////////////////

	UG_DLOG(LIB_ALG_AMG, 1, "\ngalerkin product... ");
	if(bTiming) SW.start();

	// AH = R A P
	PROFILE_BEGIN(AMGCreateAsMultiplyOf)
		CreateAsMultiplyOf(AH, R, A, PnewIndices);
	PROFILE_END();

	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

	// defragment
	if(bTiming) SW.start();
	AH.defragment();
	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, std::endl << "Defragment ... took " << SW.ms() << " ms");

#ifdef UG_PARALLEL
	AH.set_storage_type(PST_ADDITIVE);
#endif

#ifdef AMG_PRINT_AH
	UG_DLOG(LIB_ALG_AMG, 1, "AH level " << level << std::endl);
	AH.print();
#endif

	UG_DLOG(LIB_ALG_AMG, 1, "\n");

	if(m_writeMatrices)
	{
		AMGBase<TAlgebra> *t = this;
		t->write_debug_matrix_markers(level, nodes);
		this->write_debug_matrices(AH, R, AOL1, PnewIndices, level);
	}
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
#ifdef UG_PARALLEL
	m_pParallelCoarsening = NULL;
#endif
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

#ifdef UG_PARALLEL
	if(m_pParallelCoarsening)
	{	UG_LOG("Parallel Coarsening: " << m_pParallelCoarsening->tostring() << "\n"); }
	else
	{	UG_LOG("Parallel Coarsening not set.\n"); }
#endif

}

} // namespace ug

#endif //  __H__UG__LIB_ALGEBRA__RSAMG_SOLVER__RSAMG_IMPL_H__
