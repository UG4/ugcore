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
#include "amg_nodeinfo.h"

#include "amg_coarsening.h"

#include "stopwatch.h"


namespace ug{

//#define GRAPH_WITH_LOCAL_INVERSE


#define AMG_WRITE_MATRICES_PATH "/Users/mrupp/matrices/AMG_"
#define AMG_WRITE_MATRICES_MAX (200*200)

#define AMG_WRITE_GRAPH



#if 0

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
 * i is strong coupled to j, if \f$ -A_{ij} \ge theta {\max\limits}_{\phi_1(A_{ik}) < 0} { \phi_2(A_{ik}) }\f$
 * at the moment, we use the function \f$\phi_1 = \phi_2 =\f$ amg_offdiag_value
 * \param	A			matrix A for which to calculate strong connectivity graph
 * \param 	graph		the calculated strong connectivity graph of A
 *						graph is afterwards made up of connections from a node i to j if
 *						j has a strong connection to i
 * \param	theta		the smaller, the more connections are considered strong e.g. 0.25 is a reasonable default value
 * \note	you can also pre-calculate some connectivity-measurement matrix C, and call this function with C.
 */
template<typename matrix_type>
void CreateStrongConnectionGraph(const matrix_type &A, cgraph &graph, double theta)
{
	graph.resize(A.num_rows());

	for(size_t i=0; i< A.num_rows(); i++)
	{
		if(A[i].is_isolated())
			continue;

		double dmax = 0;

		for(typename matrix_type::cRowIterator conn = A.beginRow(i); !conn.isEnd(); ++conn)
		{
			if(conn.index() == i) continue; // skip diag
			if(conn.value() != 0.0 && amg_offdiag_value(conn.value()) < dmax)
				dmax = amg_offdiag_value(conn.value());
		}

		for(typename matrix_type::cRowIterator conn = A.beginRow(i); !conn.isEnd(); ++conn)
		{
			if(conn.index() == i) continue; // skip diag
			if( amg_offdiag_value(conn.value()) < theta * dmax)
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
void amg<TAlgebra>::create_AMG_level(matrix_type &AH, SparseMatrix<double> &R, const matrix_type &A,
		SparseMatrix<double> &P, int level)
{
	size_t N = A.num_rows();
	std::vector<amg_nodeinfo> nodes; nodes.resize(N);

	bool bTiming=true;
	UG_LOG("Creating level " << level << ". (" << N << " nodes)" << std::endl << std::fixed);
	stopwatch SW;

	// amg_nodeinfo: infos zu den einzelnen knoten von A, verwaltung der ratings etc

	// todo: check for isolated condition

	nodeinfo_pq_type PQ;

	//std::vector<bool coarse(N);
	int unassigned = N;

	std::vector<int> newIndex; newIndex.resize(N, -1);
	std::vector<int> posInConnections; posInConnections.resize(N, -1);


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
	graphST.create_as_transpose_of(graphS);

#ifdef AMG_WRITE_GRAPH
	WriteAMGGraphToFile(graphS, (string(AMG_WRITE_MATRICES_PATH) + "G" + ToString(level) + ".mat").c_str(), amghelper, level);
	WriteAMGGraphToFile(graphST, (string(AMG_WRITE_MATRICES_PATH) + "GT" + ToString(level) + ".mat").c_str(), amghelper, level);
#endif

	CreateMeasureOfImportancePQ(graphS, graphST, PQ, unassigned, &nodes[0]);

	if(bTiming) UG_LOG("took " << SW.ms() << " ms");

#ifdef AMG_PRINT_COARSEN_RATINGS
	 for(size_t i=0; i<N; i++)
		 UG_LOG(i << " [" << GetOriginalIndex(level, i) << "] " << nodes[i] << std::endl);
#endif

	// Coarsen
	/////////////////////////////////////////
	UG_LOG(std::endl << "coarsening... ");

	SW.start();
	int iNrOfCoarse = 0;

	Coarsen(graphST, PQ, &newIndex[0], unassigned, iNrOfCoarse, &nodes[0]);

	PreventFFConnections(graphS, graphST, &nodes[0], &newIndex[0], iNrOfCoarse);

	if(bTiming) UG_LOG("took " << SW.ms() << " ms");

	// agressive Coarsening
	/////////////////////////////////////////

	if(aggressiveCoarsening && level == 0)
	{
		// build graph 2
		//------------------

		cgraph graph2(N);

		UG_LOG(std::endl << "building graph2... ");
		if(bTiming) SW.start();

		for(size_t i=0; i < N; i++) newIndex[i] = -1;

		cgraph graphAC;
		CreateAggressiveCoarseningGraph(graphST, graphAC, &nodes[0], aggressiveCoarseningNrOfPaths, &posInConnections[0]);
		CreateMeasureOfImportanceAggressiveCoarseningPQ(graphAC, PQ, unassigned, iNrOfCoarse, &newIndex[0], &nodes[0]);
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
			Coarsen(graphAC, PQ, &newIndex[0], unassigned, iNrOfCoarse, &nodes[0]);
			//PreventFFConnections(graphS, graphST, &nodes[0], newIndex, iNrOfCoarse);
			if(bTiming) UG_LOG("took " << SW.ms() << " ms");

		}
	}

	// set parentindex for debugging
	parentIndex[level+1] = new int[iNrOfCoarse];
	for(size_t i=0; i<N; i++)
		if(nodes[i].isCoarse())
			parentIndex[level+1][ newIndex[i] ] = i;


#ifdef AMG_PRINT_COARSENING
	printCoarsening(level);
#endif

	// construct prolongation P = I_{2h->h}
	/////////////////////////////////////////

	UG_LOG(std::endl << "prolongation... ");
	unassigned = 0;

	if(bTiming) SW.start();
	CreateRugeStuebenProlongation(P, A, &newIndex[0], iNrOfCoarse, unassigned, &nodes[0], theta);
	//UG_ASSERT(unassigned == 0 || (aggressiveCoarsening && level == 0),
		//	"no aggressive coarsening, but indirect interpolation of " << unassigned " nodes ?");
	if(unassigned > 0)
		CreateIndirectProlongation(P, A, &newIndex[0], unassigned, &nodes[0], &posInConnections[0], theta);

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
	R.create_as_transpose_of(P);
	// R is already finalized

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

	// finalize
	if(bTiming) SW.start();
	AH.finalize();
	if(bTiming) UG_LOG(std::endl << "Finalizing... took " << SW.ms() << " ms");

#ifdef AMG_PRINT_AH
	UG_LOG("AH level " << level << std::endl);
	AH.print();
#endif

	UG_LOG("\n");

	delete[] posInConnections;
	delete[] newIndex;
	delete[] nodes;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// amg<matrix_type, vector_type>::init
//----------------
//! creates MG Hierachy for with matrix_type A and temporary vectors for higher levels
//! @param mat	matrix
template<typename TAlgebra>
bool amg<TAlgebra>::preprocess(matrix_type& mat)
{
	cleanup();
	A[0] = &mat;
	m_bInited = false;
	return true;
}

template<typename TAlgebra>
bool amg<TAlgebra>::init_amg()
{
	if(m_basesolver==NULL)
	{
		UG_LOG("AMG::init_amg(): No base solver selected. Call set_base_solver(basesolver) to set a base solver.\n");
		return false;
	}
	if(m_presmoother==NULL)
	{
		UG_LOG("AMG::init_amg(): No PreSmoother selected. Call set_presmoother(presmoother) to set a PreSmoother.\n");
		return false;
	}
	if(m_postsmoother==NULL)
	{
		UG_LOG("AMG::init_amg(): No PostSmoother selected. Call set_postsmoother(postsmoother) to set a PostSmoother.\n");
		return false;
	}

	// init amghelper for grid printing
	amghelper.positions = &dbg_positions[0];
	amghelper.size = A[0]->num_rows();
	amghelper.parentIndex = parentIndex;
	amghelper.dimension = dbg_dimension;

	UG_LOG("Starting AMG Setup." << std::endl << std::endl);

	stopwatch SWwhole;
	SWwhole.start();

	int i=0;
	while(i< max_levels-1)
	{
		SMO[i].setmatrix(A[i]);
		m_presmoothers[i] = m_presmoother->clone();
		m_presmoothers[i]->init(SMO[i]);
		if(m_presmoother == m_postsmoother)
			m_postsmoothers[i] = m_presmoothers[i];
		else
		{
			m_postsmoothers[i] = m_postsmoother->clone();
			m_postsmoothers[i]->init(SMO[i]);
		}

		double L = A[i]->num_rows();

#ifndef LATE_COARSE_SOLVER
		UG_LOG("No Late Coarse Solver\n");
		//if(L < 100 || A[i]->total_num_connections()/(L*L) > 0.5)	break; // abbruch falls density > 50%
		if(L < 1000)	break; // abbruch falls density > 50%
#else
		UG_LOG("Late Coarse Solver\n");
		if(L < 10)	break;
#endif
		//smoother[i].init(*A[i]);

		A[i+1] = new matrix_type;
#ifdef UG_PARALLEL
		A[i+1]->set_storage_type(PST_ADDITIVE);
#endif
		create_AMG_level(*A[i+1], R[i], *A[i], P[i], i);

		i++;
	}

	UG_ASSERT(block_traits< typename vector_type::value_type >::is_static, "dynamic not yet implemented");
	int static_nrUnknowns = block_traits< typename vector_type::value_type >::static_size;

	UG_LOG("Creating level " << i << " (" << A[i]->num_rows() << " nodes, total "
			<< A[i]->num_rows()*static_nrUnknowns << " unknowns)" << std::endl << "Using Direct Solver on Matrix "
			<< A[i]->num_rows()*static_nrUnknowns << "x" << A[i]->num_rows()*static_nrUnknowns << ". ");

	stopwatch SW; SW.start();
	SMO[i].setmatrix(A[i]);
	m_basesolver->init(SMO[i]);

	UG_LOG("Coarse Solver Setup took " << SW.ms() << "ms." << std::endl);

	used_levels = i+1;
	UG_LOG("AMG Setup finished. Used Levels: " << used_levels << ". ");
	UG_LOG("AMG Setup took " << SWwhole.ms() << " ms." << std::endl);

	// calc complexities
	double nnzs=0;
	double totallength=0;
	for(int i=0; i<used_levels; i++)
	{
		nnzs += A[i]->total_num_connections();
		totallength += A[i]->num_rows();
	}

	UG_LOG("Operator Complexity: " << nnzs/A[0]->total_num_connections() << " nodes complexity: "
			<< totallength/A[0]->num_rows() << std::endl << std::endl);

	return true;
>>>>>>> .r1641
}


//!
//! amg constructor
template<typename TAlgebra>
amg<TAlgebra>::amg() :
	amg_base<TAlgebra>(),
	eps_truncation_of_interpolation(0.3),
	theta(0.3),
	sigma(0.3),
	aggressiveCoarsening(0),
	aggressiveCoarseningNrOfPaths(2)
{
}



template<typename TAlgebra>
void amg<TAlgebra>::tostring() const
{
	amg_base<TAlgebra>::tostring();

	UG_LOG("Ruge/Stueben AMG Preconditioner:\n");
	UG_LOG("theta = " << theta << std::endl);
	UG_LOG("sigma = " << sigma << std::endl);

	if(aggressiveCoarsening)	{UG_LOG("Aggressive Coarsening is on, A" << aggressiveCoarseningNrOfPaths << "-mode." << std::endl);}
	else						{UG_LOG("no Aggressive Coarsening" << std::endl);}
}

} // namespace ug

#endif //  __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_IMPL_H__
