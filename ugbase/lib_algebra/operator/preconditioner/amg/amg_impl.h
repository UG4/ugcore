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
 * \param	m_dTheta		the smaller, the more connections are considered strong e.g. 0.25 is a reasonable default value
 * \note	you can also pre-calculate some connectivity-measurement matrix C, and call this function with C.
 */
template<typename matrix_type>
void CreateStrongConnectionGraph(const matrix_type &A, cgraph &graph, double m_dTheta)
{
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
	size_t N = A.num_rows();
	stdvector<amg_nodeinfo> nodes; nodes.resize(N);

	bool bTiming=true;
	UG_LOG("Creating level " << level << ". (" << N << " nodes)" << std::endl << std::fixed);
	stopwatch SW;

	// amg_nodeinfo: infos zu den einzelnen knoten von A, verwaltung der ratings etc

	// todo: check for isolated condition

	nodeinfo_pq_type PQ;

	//stdvector<bool coarse(N);
	int unassigned = N;

	stdvector<int> newIndex; newIndex.resize(N, -1);
	stdvector<int> posInConnections; posInConnections.resize(N, -1);


	UG_LOG("A.totalNrOfConnections = " << A.total_num_connections() << std::endl);

	cgraph graphS(N);
	cgraph graphST;

	// build graph
	/////////////////////////////////////////

	UG_LOG("building graph... ");
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

	CreateMeasureOfImportancePQ(graphS, graphST, PQ, unassigned, nodes);

	if(bTiming) UG_LOG("took " << SW.ms() << " ms");

	IF_DEBUG(LIB_ALG_AMG, 4)
	{
		UG_LOG("Coarsen ratings:\n")
		for(size_t i=0; i<N; i++)
			UG_LOG(i << " [" << m_amghelper.GetOriginalIndex(level, i) << "] " << nodes[i] << std::endl);
	}

	// Coarsen
	/////////////////////////////////////////
	UG_LOG(std::endl << "coarsening... ");

	SW.start();
	int iNrOfCoarse = 0;

	Coarsen(graphST, PQ, newIndex, unassigned, iNrOfCoarse, nodes);

	PreventFFConnections(graphS, graphST, nodes, newIndex, iNrOfCoarse);

	if(bTiming) UG_LOG("took " << SW.ms() << " ms");

	// agressive Coarsening
	/////////////////////////////////////////

	if(m_bAggressiveCoarsening && level == 0)
	{
		// build graph 2
		//------------------

			UG_LOG(std::endl << "building graph2... ");
		if(bTiming) SW.start();

		for(size_t i=0; i < N; i++) newIndex[i] = -1;

		unassigned = 0;
		cgraph graphAC(N);
		CreateAggressiveCoarseningGraph(graphST, graphAC, nodes, m_iAggressiveCoarseningNrOfPaths, &posInConnections[0]);
		CreateMeasureOfImportanceAggressiveCoarseningPQ(graphAC, PQ, unassigned, iNrOfCoarse, newIndex, nodes);
		if(bTiming) UG_LOG("took " << SW.ms() << " ms");
		// coarsen 2
		//------------------

		if(unassigned == 0)
		{
			UG_LOG(std::endl << "skipping coarsening2: no unassigned nodes.");
		}
		else
		{
			UG_LOG(std::endl << "coarsening2... ");

			if(bTiming) SW.start();
			Coarsen(graphAC, PQ, newIndex, unassigned, iNrOfCoarse, nodes);
			//PreventFFConnections(graphS, graphST, nodes, newIndex, iNrOfCoarse);
			if(bTiming) UG_LOG("took " << SW.ms() << " ms");

		}
	}

	// set parentindex for debugging
	m_parentIndex.resize(level+2);
	m_parentIndex[level+1].resize(iNrOfCoarse);
	for(size_t i=0; i<N; i++)
		if(nodes[i].isCoarse())
			m_parentIndex[level+1][ newIndex[i] ] = i;


#ifdef AMG_PRINT_COARSENING
	printCoarsening(level);
#endif

	// construct prolongation P = I_{2h->h}
	/////////////////////////////////////////

	UG_LOG(std::endl << "prolongation... ");
	unassigned = 0;

	if(bTiming) SW.start();
	CreateRugeStuebenProlongation(P, A, newIndex, iNrOfCoarse, unassigned, nodes, m_dTheta);
	//UG_ASSERT(unassigned == 0 || (m_bAggressiveCoarsening && level == 0),
		//	"no aggressive coarsening, but indirect interpolation of " << unassigned " nodes ?");
	if(unassigned > 0)
		CreateIndirectProlongation(P, A, newIndex, unassigned, nodes, &posInConnections[0], m_dTheta);

	P.defragment();

/*#ifdef UG_PARALLEL
	P.set_storage_type(PST_CONSISTENT);
#endif*/

	if(bTiming) UG_LOG("took " << SW.ms() << " ms");

#ifdef AMG_PRINT_P
	UG_LOG(std::endl << "Prolongation level " << level << std::endl);
	P.print();
#endif


	// construct restriction R = I_{h->2h}
	/////////////////////////////////////////

	UG_LOG(std::endl << "restriction... ");
	if(bTiming) SW.start();

	// construct restriction R = I_{h -> 2h}
	R.set_as_transpose_of(P);
	// R is already defragmented

/*#ifdef UG_PARALLEL
	R.set_storage_type(PST_CONSISTENT);
#endif*/

	if(bTiming) UG_LOG("took " << SW.ms() << " ms");

#ifdef AMG_PRINT_R
	UG_LOG(std::endl << "Restriction level " << level << std::endl);
	R.print();
#endif

	// create Galerkin product
	/////////////////////////////////////////

	UG_LOG("\ngalerkin product... ");
	if(bTiming) SW.start();

	// AH = R A P
	CreateAsMultiplyOf(AH, R, A, P);

	if(bTiming) UG_LOG("took " << SW.ms() << " ms");

	// defragment
	if(bTiming) SW.start();
	AH.defragment();
	if(bTiming) UG_LOG(std::endl << "Defragment ... took " << SW.ms() << " ms");

#ifdef AMG_PRINT_AH
	UG_LOG("AH level " << level << std::endl);
	AH.print();
#endif

	UG_LOG("\n");


	if(m_writeMatrices && A.num_rows() < AMG_WRITE_MATRICES_MAX)
	{
		std::fstream ffine((std::string(m_writeMatrixPath) + "AMG_fine" + ToString(level) + ".marks").c_str(), std::ios::out);
		std::fstream fcoarse((std::string(m_writeMatrixPath) + "AMG_coarse" + ToString(level) + ".marks").c_str(), std::ios::out);
		std::fstream fother((std::string(m_writeMatrixPath) + "AMG_other" + ToString(level) + ".marks").c_str(), std::ios::out);
		std::fstream fdirichlet((std::string(m_writeMatrixPath) + "AMG_dirichlet" + ToString(level) + ".marks").c_str(), std::ios::out);
		for(size_t i=0; i<N; i++)
		{
			int o = m_amghelper.GetOriginalIndex(level, i);
			if(nodes[i].isFineDirect()) ffine << o << "\n";
			else if(nodes[i].isCoarse()) fcoarse << o << "\n";
			else fother << o << "\n";
		}

		UG_LOG("write matrices");
		AMGWriteToFile(A, level, level, (std::string(m_writeMatrixPath) + "AMG_A" + ToString(level) + ".mat").c_str(), m_amghelper);
		std::fstream f((std::string(m_writeMatrixPath) + "AMG_A" + ToString(level) + ".mat").c_str(), std::ios::out | std::ios::app);
		f << "c " << std::string(m_writeMatrixPath) << "AMG_fine" << level << ".marks\n";
		f << "c " << std::string(m_writeMatrixPath) << "AMG_coarse" << level << ".marks\n";
		f << "c " << std::string(m_writeMatrixPath) << "AMG_other" << level << ".marks\n";
		f << "c " << std::string(m_writeMatrixPath) << "AMG_dirichlet" << level << ".marks\n";
		f << "v " << std::string(m_writeMatrixPath) << "AMG_d" << level << ".values\n";
		UG_LOG(".");

		AMGWriteToFile(P, level+1, level, (std::string(m_writeMatrixPath) + "AMG_Pp" + ToString(level) + ".mat").c_str(), m_amghelper);
		std::fstream f2((std::string(m_writeMatrixPath) + "AMG_Pp" + ToString(level) + ".mat").c_str(), std::ios::out | std::ios::app);
		f2 << "c " << std::string(m_writeMatrixPath) << "AMG_fine" << level << ".marks\n";
		f2 << "c " << std::string(m_writeMatrixPath) << "AMG_coarse" << level << ".marks\n";
		f2 << "c " << std::string(m_writeMatrixPath) << "AMG_other" << level << ".marks\n";
		f2 << "c " << std::string(m_writeMatrixPath) << "AMG_dirichlet" << level << ".marks\n";
		UG_LOG(".");

		AMGWriteToFile(R, level, level+1, (std::string(m_writeMatrixPath) + "AMG_Rr" + ToString(level) + ".mat").c_str(), m_amghelper);
		std::fstream f3((std::string(m_writeMatrixPath) + "AMG_Rr" + ToString(level) + ".mat").c_str(), std::ios::out | std::ios::app);
		f3 << "c " << std::string(m_writeMatrixPath) << "AMG_fine" << level << ".marks\n";
		f3 << "c " << std::string(m_writeMatrixPath) << "AMG_coarse" << level << ".marks\n";
		f3 << "c " << std::string(m_writeMatrixPath) << "AMG_other" << level << ".marks\n";
		f3 << "c " << std::string(m_writeMatrixPath) << "AMG_dirichlet" << level << ".marks\n";
		UG_LOG(".");

		AMGWriteToFile(AH, level+1, level+1, (std::string(m_writeMatrixPath) + "AMG_A" + ToString(level+1) + ".mat").c_str(), m_amghelper);
		UG_LOG(". done.\n");
	}
}


//!
//! amg constructor
template<typename TAlgebra>
amg<TAlgebra>::amg() :
	amg_base<TAlgebra>(),
	m_dEpsilon(0.3),
	m_dTheta(0.3),
	m_dSigma(0.3),
	m_bAggressiveCoarsening(0),
	m_iAggressiveCoarseningNrOfPaths(2)
{
}



template<typename TAlgebra>
void amg<TAlgebra>::tostring() const
{
	amg_base<TAlgebra>::tostring();

	UG_LOG("Ruge/Stueben AMG Preconditioner:\n");
	UG_LOG(" Theta (strong connectivity) = " << m_dTheta << std::endl);
	UG_LOG(" Sigma = " << m_dSigma << std::endl);
	UG_LOG(" Epsilon (truncation of interpolation) = " << m_dEpsilon << std::endl);

	if(m_bAggressiveCoarsening)	{UG_LOG(" Aggressive Coarsening is on, A" << m_iAggressiveCoarseningNrOfPaths << "-mode." << std::endl);}
	else						{UG_LOG(" no Aggressive Coarsening" << std::endl);}
}

} // namespace ug

#endif //  __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_IMPL_H__
