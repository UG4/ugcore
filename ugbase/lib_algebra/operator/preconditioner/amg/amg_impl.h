/**
 * \file amg_impl.h
 *
 * \author Martin Rupp
 *
 * \date 03.12.2009
 *
 * implementation file for amg
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#ifndef __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_IMPL_H__
#define __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_IMPL_H__

//#include "sparsematrix_util.h"

#include "amg.h"
#include "graph.h"

#include "amg_rs_prolongation.h"
#include "amg_debug.h"
#include "amg_debug_helper.h"
#include "amg_nodeinfo.h"

#include "amg_coarsening.h"

#include "stopwatch.h"
#include "amg_profiling.h"

#ifdef UG_PARALLEL
#include "pcl/pcl.h"
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
void amg<TAlgebra>::FullSubdomainBlocking(matrix_type &AH, IndexLayout &nextMasterLevelLayout, IndexLayout &nextSlaveMasterLayout)
{

}*/

template<typename TAlgebra>
void amg<TAlgebra>::create_new_indices(stdvector<int> &newIndex, const AMGNodes &nodes, size_t level)
{
	AMG_PROFILE_FUNC();
	is_fine.resize(level+1);
	is_fine[level].resize(nodes.size());
	size_t iNrOfCoarse=0;
	for(size_t i=0; i < nodes.size(); i++)
		if(nodes[i].is_coarse())
		{
			newIndex[i] = iNrOfCoarse;
			is_fine[level][i] = false;
			iNrOfCoarse++;
		}
		else
		{
			is_fine[level][i] = true;
			newIndex[i] = -1;
		}
	UG_ASSERT(iNrOfCoarse == nodes.get_nr_of_coarse(), "");
}

template<typename TAlgebra>
void amg<TAlgebra>::create_parentIndex(const stdvector<int> &newIndex, const AMGNodes &nodes, size_t level)
{
	AMG_PROFILE_FUNC();
	m_parentIndex.resize(level+2);
	m_parentIndex[level+1].resize(nodes.get_nr_of_coarse());
	for(size_t i=0; i<nodes.size(); i++)
		if(nodes[i].is_coarse())
			m_parentIndex[level+1][ newIndex[i] ] = i;
}



template<typename TAlgebra>
void amg<TAlgebra>::debug_matrix_write(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
		prolongation_matrix_type &P, size_t level, const AMGNodes &nodes)
{
	AMG_PROFILE_FUNC();
	std::fstream ffine((std::string(m_writeMatrixPath) + "AMG_fine_L" + ToString(level) + ".marks").c_str(), std::ios::out);
	std::fstream fcoarse((std::string(m_writeMatrixPath) + "AMG_coarse_L" + ToString(level) + ".marks").c_str(), std::ios::out);
	std::fstream fother((std::string(m_writeMatrixPath) + "AMG_other_L" + ToString(level) + ".marks").c_str(), std::ios::out);
	std::fstream fdirichlet((std::string(m_writeMatrixPath) + "AMG_dirichlet_L" + ToString(level) + ".marks").c_str(), std::ios::out);
	for(size_t i=0; i<nodes.size(); i++)
	{
		int o = m_amghelper.GetOriginalIndex(level, i);
		if(nodes[i].is_fine_direct()) ffine << o << "\n";
		else if(nodes[i].is_coarse()) fcoarse << o << "\n";
		else fother << o << "\n";
	}

	UG_DLOG(LIB_ALG_AMG, 1, "write matrices");
	AMGWriteToFile(A, level, level, (std::string(m_writeMatrixPath) + "AMG_A_L" + ToString(level) + ".mat").c_str(), m_amghelper);
	std::fstream f((std::string(m_writeMatrixPath) + "AMG_A_L" + ToString(level) + ".mat").c_str(), std::ios::out | std::ios::app);
	f << "c " << std::string(m_writeMatrixPath) << "AMG_fine_L" << level << ".marks\n";
	f << "c " << std::string(m_writeMatrixPath) << "AMG_coarse_L" << level << ".marks\n";
	f << "c " << std::string(m_writeMatrixPath) << "AMG_other_L" << level << ".marks\n";
	f << "c " << std::string(m_writeMatrixPath) << "AMG_dirichlet_L" << level << ".marks\n";
	f << "v " << std::string(m_writeMatrixPath) << "AMG_d_L" << level << ".values\n";
	UG_DLOG(LIB_ALG_AMG, 1, ".");

	AMGWriteToFile(P, level+1, level, (std::string(m_writeMatrixPath) + "AMG_Pp_L" + ToString(level) + ".mat").c_str(), m_amghelper);
	std::fstream f2((std::string(m_writeMatrixPath) + "AMG_Pp_L" + ToString(level) + ".mat").c_str(), std::ios::out | std::ios::app);
	f2 << "c " << std::string(m_writeMatrixPath) << "AMG_fine_L" << level << ".marks\n";
	f2 << "c " << std::string(m_writeMatrixPath) << "AMG_coarse_L" << level << ".marks\n";
	f2 << "c " << std::string(m_writeMatrixPath) << "AMG_other_L" << level << ".marks\n";
	f2 << "c " << std::string(m_writeMatrixPath) << "AMG_dirichlet_L" << level << ".marks\n";
	UG_DLOG(LIB_ALG_AMG, 1, ".");

	AMGWriteToFile(R, level, level+1, (std::string(m_writeMatrixPath) + "AMG_Rr_L" + ToString(level) + ".mat").c_str(), m_amghelper);
	std::fstream f3((std::string(m_writeMatrixPath) + "AMG_Rr_L" + ToString(level) + ".mat").c_str(), std::ios::out | std::ios::app);
	f3 << "c " << std::string(m_writeMatrixPath) << "AMG_fine_L" << level << ".marks\n";
	f3 << "c " << std::string(m_writeMatrixPath) << "AMG_coarse_L" << level << ".marks\n";
	f3 << "c " << std::string(m_writeMatrixPath) << "AMG_other_L" << level << ".marks\n";
	f3 << "c " << std::string(m_writeMatrixPath) << "AMG_dirichlet_L" << level << ".marks\n";
	UG_DLOG(LIB_ALG_AMG, 1, ".");

	AMGWriteToFile(AH, level+1, level+1, (std::string(m_writeMatrixPath) + "AMG_A_L" + ToString(level+1) + ".mat").c_str(), m_amghelper);
	UG_DLOG(LIB_ALG_AMG, 1, ". done.\n");
}

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
void amg<TAlgebra>::create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
		prolongation_matrix_type &P, size_t level)
{
	AMG_PROFILE_FUNC();
#ifdef UG_PARALLEL
	UG_ASSERT(pcl::GetNumProcesses()==1, "not implemented for procs > 1");
#endif
	size_t N = A.num_rows();
	AMGNodes nodes(N);

	bool bTiming=true;
	UG_DLOG(LIB_ALG_AMG, 1, "Creating level " << level << ". (" << N << " nodes)" << std::endl << std::fixed);
	stopwatch SW;

	// todo: check for isolated condition

	nodeinfo_pq_type PQ;

	stdvector<int> posInConnections; posInConnections.resize(N, -1);

	UG_DLOG(LIB_ALG_AMG, 1, "A.totalNrOfConnections = " << A.total_num_connections() << std::endl);

	cgraph graphS(N);
	cgraph graphST;

	// build graph
	/////////////////////////////////////////

	UG_DLOG(LIB_ALG_AMG, 1, "building graph... ");
	if(bTiming) SW.start();

	CreateStrongConnectionGraph(A, graphS);

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

	Coarsen(graphST, PQ, nodes);

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

	// create new indices
	stdvector<int> newIndex(nodes.size());
	create_new_indices(newIndex, nodes, level);


	// create parentIndex (for debug only)
	if(m_writeMatrices)
		create_parentIndex(newIndex, nodes, level);

	// construct prolongation P = I_{2h->h}
	/////////////////////////////////////////

	UG_DLOG(LIB_ALG_AMG, 1, std::endl << "prolongation... ");


	if(bTiming) SW.start();
	CreateRugeStuebenProlongation(P, A, nodes, newIndex, m_dTheta);
	UG_ASSERT(nodes.get_unassigned() == 0 || (m_bAggressiveCoarsening && level == 0),
			"no aggressive coarsening, but indirect interpolation of " << nodes.get_unassigned() << " nodes ?");
	if(nodes.get_unassigned() > 0)
		CreateIndirectProlongation(P, A, nodes, &posInConnections[0], m_dTheta);

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
#endif

#ifdef AMG_PRINT_AH
	UG_DLOG(LIB_ALG_AMG, 1, "AH level " << level << std::endl);
	AH.print();
#endif

	UG_DLOG(LIB_ALG_AMG, 1, "\n");

	if(m_writeMatrices && A.num_rows() < AMG_WRITE_MATRICES_MAX)
		debug_matrix_write(AH, R, A, P, level, nodes);
}


//!
//! amg constructor
template<typename TAlgebra>
amg<TAlgebra>::amg() :
	amg_base<TAlgebra>(),
	m_dEpsilonTr(0.3),
	m_dTheta(0.3),
	//m_dSigma(0.3),
	m_bAggressiveCoarsening(0),
	m_iAggressiveCoarseningNrOfPaths(2)
{
}



template<typename TAlgebra>
void amg<TAlgebra>::tostring() const
{
	amg_base<TAlgebra>::tostring();

	UG_LOG("Ruge/Stueben AMG Preconditioner:\n");
	UG_LOG(" theta (epsilon_strong, strong connectivity) = " << m_dTheta << std::endl);
	//UG_LOG(" sigma = " << m_dSigma << std::endl);
	UG_LOG(" epsilon_tr (truncation of interpolation) = " << m_dEpsilonTr << std::endl);

	if(m_bAggressiveCoarsening)	{UG_LOG(" Aggressive Coarsening is on, A" << m_iAggressiveCoarseningNrOfPaths << "-mode." << std::endl);}
	else						{UG_LOG(" no Aggressive Coarsening" << std::endl);}
}

} // namespace ug

#endif //  __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_IMPL_H__
