/**
 * \file famg.cpp
 *
 * \author Martin Rupp
 *
 * \date 7.12.2010
 *
 * implementation file for famg
 *
 * Goethe-Center for Scientific Computing 2010-2011.
 *
 * NOTE: for test purposes, functions are here in a cpp file.
 */


#ifndef __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__
#define __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__

//#include "sparsematrix_util.h"

#include "common/assert.h"
#include "common/log.h"

#include "famg.h"
#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/cpu_algebra_types.h"
#include "lib_algebra/common/stl_debug.h"
#include "../amg_debug_helper.h"

#include "../stopwatch.h"
#include "common/assert.h"
#include "../maxheap.h"
#include "famg.h"
#include <set>
#include <string>
#include <stack>

// for external coarsening
#include "../rsamg/rsamg.h"
#include "../rsamg/rsamg_coarsening.h"

#ifdef UG_PARALLEL
#include "lib_algebra/parallelization/parallel_matrix_overlap_impl.h"
#endif

#include "lib_algebra/common/connection_viewer_output.h"

std::stack<int> g_DebugLevelStack;



#ifdef UG_DEBUG
#define GET_DEBUG_LEVEL(tag) ug::GetLogAssistant().get_debug_level(ug::LogAssistant::tag)
#else
#define GET_DEBUG_LEVEL(tag) 0
#endif
#define PUSH_DEBUG_LEVEL() g_DebugLevelStack.push(GET_DEBUG_LEVEL(LIB_ALG_AMG))
#define SET_AND_PUSH_DEBUG_LEVEL(level) g_DebugLevelStack.push(GET_DEBUG_LEVEL(LIB_ALG_AMG)); ug::GetLogAssistant().set_debug_level(level)
#define POP_DEBUG_LEVEL() ug::GetLogAssistant().set_debug_level(g_DebugLevelStack.pop())


#ifdef UG_PARALLEL
std::string GetProcFilename(std::string name, std::string extension)
{
	return name + "_" + ToString(pcl::GetProcRank()) + extension;
}

#else
std::string GetProcFilename(std::string name, std::string extension)
{
	return name + "_0" + extension;
}
#endif
std::string GetProcFilename(std::string path, std::string name, std::string extension)
{
	return path + GetProcFilename(name, extension);
}




template<typename vector_type>
void print_vector(const vector_type &vec, const char *p)
{
	UG_LOG(p << ": ");
	for(size_t i=0; i<vec.size(); i++)
	{
		if(i>0) UG_LOG(", ");
		UG_LOG(vec[i]);
	}
	UG_LOG("\n");
}

//const double damping = 0.5;

#if 1

//#define AMG_WRITE_GRAPH


//#define AMG_PRINT_INDIRECT

#define FAMG_PRINT_GRAPH

#define FAMG_PRINT_COARSENING
#define FAMG_PRINT_P
#define FAMG_PRINT_R
#define FAMG_PRINT_AH

#define FAMG_PRINT_POSSIBLE_PARENTS
#endif

//#define GRAPH_WITH_LOCAL_INVERSE
#include "famg_structs.h"
#include "famg_interpolation_calculator.h"
#include "famg_nodeinfo.h"

///	algebra blocks are static.
template <> struct block_traits<MathVector<3> >
{
	enum{
		is_static = 1
	};
};


// CreateSymmConnectivityGraph:
//------------------------------
/**
 * Creates a graph which has a connection from i to j if A_{ij} != 0 or A_{ji} != 0
 * \param A						the matrix A
 * \param SymmNeighGraph		later used to determine which neighbors' rating needs to be updated
 */
template<typename matrix_type>
void CreateSymmConnectivityGraph(const matrix_type &A, cgraph &SymmNeighGraph)
{
	AMG_PROFILE_FUNC();

	for(size_t i=0; i<A.num_rows(); i++)
		for(typename matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
		{
			if(i == conn.index() || conn.value() == 0.0) continue;
			SymmNeighGraph.set_connection(i, conn.index());
			SymmNeighGraph.set_connection(conn.index(), i);
		}

#ifdef FAMG_PRINT_GRAPH

	IF_DEBUG(LIB_ALG_AMG, 3)
	{
		UG_LOG("\nSymmNeighGraph:\n");
		SymmNeighGraph.print();
		UG_LOG("\n\n");
	}
#endif
}

template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
class FAMGLevelCalculator
{
private:
	// refernce to famg object. this will be done differently in the future
	// (constructor looks ugly)
	FAMG<CPUAlgebra> &m_famg;

	// the coarse matrices to be filled
	matrix_type &AH;

	// the fine grid matrix (const!)
	const matrix_type &A;

	bool bTiming; //< if true, timing will be outputed

	// interpolation/restriction matrices to be filled
	prolongation_matrix_type &R, &PnewIndices;

	size_t level; // current amg level

	// the calculator which calculates our interpolation weighs
	FAMGInterpolationCalculator<matrix_type, vector_type> calculator;

	// this structure holds coarse/fine information as well as
	// master/slave. it tries to replace the missing "node" object here
	FAMGNodes rating;

	//! list of possible interpolating parents for each node
	stdvector<stdvector<neighborstruct2> > possible_parents;

	stdvector<bool> prolongation_calculated;

	//! heap used for sorting of ratings
	maxheap<FAMGNode> heap;

	//! used to determine which neighbors' rating needs to be updated:
	cgraph SymmNeighGraph;

	// our testvectors
	stdvector< vector_type > &m_testvectors;

#ifdef UG_PARALLEL
	//! layout for sending/receiving coarsening information
	//! note: this is not a classical master/slave distribution!
	IndexLayout OLCoarseningReceiveLayout, OLCoarseningSendLayout;

	//! processes with higher/lower color to distribute coarsening information
	std::vector<int> processesWithLowerColor, processesWithHigherColor;
	int m_myColor;
#endif

#ifdef UG_PARALLEL
	// overlap 2 matrix in parallel case
	matrix_type A_OL2;
	IndexLayout &nextLevelMasterLayout;
	IndexLayout &nextLevelSlaveLayout;

#else
	const matrix_type &A_OL2;
#endif

	prolongation_matrix_type PoldIndices;


	void rs_amg_external_coarsening()
	{
		AMG_PROFILE_FUNC();
		stopwatch SW;
		UG_DLOG(LIB_ALG_AMG, 1, std::endl << "external coarsening... "); if(bTiming) SW.start();
		nodeinfo_pq_type PQ;
		size_t N = A.num_rows();

		// use RSAMG's methods to do standard rs coarsening
		AMGNodes nodes(N);
		cgraph graphS, graphST;
		CreateStrongConnectionGraph(A, graphS, 0.3);

		graphST.set_as_transpose_of(graphS);
		CreateMeasureOfImportancePQ(graphS, graphST, PQ, nodes);
		Coarsen(graphST, PQ, nodes);
		PreventFFConnections(graphS, graphST, nodes);


		if(m_famg.m_writeMatrices)
			write_debug_matrix_markers();

		// aggressive coarsening
		if(0 && m_famg.get_aggressive_coarsening() == true && level == 0)
		{
			UG_DLOG(LIB_ALG_AMG, 2, std::endl << "building graph2... ");
			UG_DLOG(LIB_ALG_AMG, 2, "unassigned = " << nodes.get_unassigned() << "\n");

			//unassigned = 0; ??
			cgraph graphAC(N);
			size_t m_iAggressiveCoarseningNrOfPaths = 2;
			stdvector<int> posInConnections; posInConnections.resize(N, -1);
			CreateAggressiveCoarseningGraph(graphST, graphAC, nodes, m_iAggressiveCoarseningNrOfPaths, &posInConnections[0]);
			CreateMeasureOfImportanceAggressiveCoarseningPQ(graphAC, PQ, nodes);
			if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");
			// coarsen 2
			//------------------

			if(nodes.get_unassigned() == 0)
			{
				UG_DLOG(LIB_ALG_AMG, 2, std::endl << "skipping coarsening2: no unassigned nodes.");
			}
			else
			{
				UG_DLOG(LIB_ALG_AMG, 2, std::endl << "coarsening2... ");
				Coarsen(graphAC, PQ, nodes);
				//PreventFFConnections(graphS, graphST, nodes);
			}
		}

		// coarsening done. transfer AMG nodeinfo -> FAMG nodeinfo
		for(size_t i=0; i<N; i++)
		{
			if(nodes[i].is_coarse())
				rating.set_coarse(i);
		}

		for(size_t i=0; i<N; i++)
		{
			if(nodes[i].is_fine_direct())
			{
				if(graphS.is_isolated(i))
					rating.set_fine(i);
				else
					calculator.get_all_neighbors_interpolation(i, PoldIndices, rating);
			}
		}

		if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");
		//if(nodes[i].is_unassigned_fine_indirect())
		//	rating.set_aggressive_fine(i);
	}

	void create_restriction()
	{
		AMG_PROFILE_FUNC();
		stopwatch SW;
		UG_DLOG(LIB_ALG_AMG, 1, std::endl << "restriction... "); if(bTiming) SW.start();
		// construct restriction R = I_{h -> 2h}
		R.set_as_transpose_of(PnewIndices);
		// R is already defragmented
		if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms.");

	#ifdef UG_PARALLEL
		R.set_storage_type(PST_CONSISTENT);
	#endif

	#ifdef FAMG_PRINT_R
		IF_DEBUG(LIB_ALG_AMG, 4)
		{
			UG_DLOG(LIB_ALG_AMG, 1, std::endl << "Restriction level " << level << std::endl);
			R.print();
		}
	#endif
	}

	void create_galerkin_product()
	{
		AMG_PROFILE_FUNC();
		stopwatch SW;
		UG_DLOG(LIB_ALG_AMG, 1, "\ngalerkin product... "); if(bTiming) SW.start();
		// AH = R A P

		/* matrix_type AH2;
		matrix_type RoldIndices;
		RoldIndices.set_as_transpose_of(PoldIndices);
		CreateAsMultiplyOf2(AH2, RoldIndices, A, PoldIndices); */

		// R.print("R");
		// A.print("A");
		// PnewIndices.print("P");
		CreateAsMultiplyOf(AH, R, A, PnewIndices, 1e-12);
		// AH.print();
		if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

		// finalize
		if(bTiming) { UG_DLOG(LIB_ALG_AMG, 1, std::endl << "Finalizing.."); SW.start(); }
		AH.defragment();
		if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, " took " << SW.ms() << " ms\n");



	#ifdef UG_PARALLEL
		AH.set_storage_type(PST_ADDITIVE);
		AH.set_communicator(A_OL2.get_communicator());
		AH.set_process_communicator(A_OL2.get_process_communicator());
		// AH.set_layouts(nextLevelMasterLayout, nextLevelSlaveLayout); // is already set
	#endif


		PROFILE_END();

	#ifdef FAMG_PRINT_AH
		IF_DEBUG(LIB_ALG_AMG, 3)
		{
			UG_DLOG(LIB_ALG_AMG, 1, "AH level " << level << std::endl);
			AH.print();
		}
	#endif
	}
public:
	// todo: clean up this mess of constructor
	FAMGLevelCalculator(FAMG<CPUAlgebra> &f,
			matrix_type &_AH, prolongation_matrix_type &_R,  const matrix_type &_A,
			prolongation_matrix_type &_P, size_t _level,
			stdvector< vector_type > &testvectors, stdvector<double> &omega)
	: m_famg(f), AH(_AH), A(_A), R(_R), PnewIndices(_P), level(_level),
				calculator(A, A_OL2,
					m_famg.get_delta(),
					m_famg.get_theta(),
					m_famg.get_damping_for_smoother_in_interpolation_calculation(),
					m_famg.get_epsilon_truncation(),
					testvectors, omega),
			rating(PoldIndices, _level, m_famg.m_amghelper),
			m_testvectors(testvectors)
#ifndef UG_PARALLEL
			, A_OL2(A)
#else
			, nextLevelMasterLayout(AH.get_master_layout()), nextLevelSlaveLayout(AH.get_slave_layout())
#endif
	{
	}

public:
	void do_calculation()
	{
		AMG_PROFILE_FUNC();

		UG_DLOG(LIB_ALG_AMG, 1, "Creating level " << level << ". (" << A.num_rows() << " nodes)" << std::fixed);
		bTiming=true;
		stopwatch SW, SWwhole; SWwhole.start();

		// 1. Overlap calculation
		//---------------------------
#ifdef UG_PARALLEL
		create_OL2_matrix();
#else
		rating.create(A_OL2.num_rows());
#endif
		size_t N = A_OL2.num_rows();

		// 2. global Testvector calculation (damping)
		//-----------------------------------------------

		calculate_testvectors();
		calculator.init();

		// 3. create heap, P, SymmNeighGraph
		//-------------------------------------

		bool bExternalCoarsening = m_famg.m_bExternalCoarsening;
		bool bUsePrecalculate = m_famg.m_bUsePrecalculate;
		if(bExternalCoarsening == false)
		{
			UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, m_famg.iDebugLevelPhase3);
			heap.create(rating.nodes);

			IF_DEBUG(LIB_ALG_AMG, 3)
			{
				UG_DLOG(LIB_ALG_AMG, 1, "\bThe Matrix A:\n\n");
				A.print();
				UG_DLOG(LIB_ALG_AMG, 1, "\n\n");
			}

			UG_DLOG(LIB_ALG_AMG, 1, "\nCreate P, SymmNeighGraph... "); if(bTiming) SW.start();
			PoldIndices.resize(N, N);
			// get neighboring information
			SymmNeighGraph.resize(N);
			CreateSymmConnectivityGraph(A_OL2, SymmNeighGraph);
			if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

			if(bUsePrecalculate)
				// get possible parent nodes
				calculate_all_possible_parent_pairs();
		}
		else
			PoldIndices.resize(N, N);


	#ifdef UG_PARALLEL
		// 4. coloring in parallel
		//-------------------------------------
		color_process_graph();
		receive_coarsening_from_processes_with_lower_color();
	#endif

		rating.calculate_unassigned();

		if(bExternalCoarsening == false)
		{
			// 5. do coarsening
			//-------------------------
			if(bUsePrecalculate)
			{
				// calculate ratings (not precalculateable because of coarse/uninterpolateable)
				UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, m_famg.iDebugLevelGetRatings);
				UG_DLOG(LIB_ALG_AMG, 1, std::endl << "calculate ratings... "); if(bTiming) SW.start();
				GetRatings(possible_parents, rating, heap);
				if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

				//heap.print();

				// do coarsening
				precalculate_coarsening();
			}
			else
				on_demand_coarsening();


			// 7. do aggressive coarsening
			//-----------------------------------------------
			// this has to be BEFORE communication of coarse/fine markers
			//

			UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, m_famg.iDebugLevelAggressiveCoarsening);
			UG_DLOG(LIB_ALG_AMG, 1, std::endl << "second coarsening... "); if(bTiming) SW.start();
			if(m_famg.get_aggressive_coarsening() == true)
				get_aggressive_coarsening_interpolation();
			else
				set_uninterpolateable_as_coarse();
			if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms.");

			IF_DEBUG(LIB_ALG_AMG, 3)
			{	rating.print(); 	}

			#ifdef UG_PARALLEL
			// 6. send coarsening data to processes with higher color
			//----------------------------------------------------------
				send_coarsening_data_to_processes_with_higher_color();
				// we could use this information to do better interpolation on the border
				// send_coarsening_data_to_processes_with_lower_color_and_receive_coarsening_data_from_processes_with_higher_color

				// 7. calculate interpolation of nodes with uncalculated_fine
				//calculate_uncalculated_fine_nodes();
			#endif

			// 8.

			// resize P, output statistics & debug stuff

			UG_DLOG(LIB_ALG_AMG, 1, std::endl << N - rating.get_unassigned() << " nodes assigned, " << rating.get_nr_of_coarse() << " coarse, "
					<< N - rating.get_unassigned() - rating.get_nr_of_coarse() << " fine, " << rating.get_unassigned() << " unassigned.");

			// set P to real size

			#ifdef UG_PARALLEL
				communicate_prolongation();
			#endif

			// [ debug output
			if(m_famg.m_writeMatrices)
				write_debug_matrix_markers();
			// ]
		}
		else
			rs_amg_external_coarsening();

		UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, m_famg.iDebugLevelAfterCommunicateProlongation);
		IF_DEBUG(LIB_ALG_AMG, 3)
		{	rating.print(); 	}

		create_new_indices();
		create_fine_marks();

		//UG_DLOG(LIB_ALG_AMG, 1, "parentIndex level " << level << "\n")
		//for(size_t i=0; i<rating.get_nr_of_coarse(); i++) { UG_ASSERT(amghelper.parentIndex[level+1][i] != -1, i << " == -1???"); UG_LOG(i << " = " << amghelper.parentIndex[level+1][i] << "\n"); }

		UG_DLOG(LIB_ALG_AMG, 1, std::endl << rating.get_nr_of_coarse() << " coarse, " << N - rating.get_unassigned() - rating.get_nr_of_coarse() << " fine.");

		// 8. construct restriction R = I_{h->2h}
		//-----------------------------------------
		create_restriction();

		// 9. create Galerkin product AH = R A P
		//-----------------------------------------

		create_galerkin_product();
		// 10.
		//-----------------------------------------

		if(m_famg.m_writeMatrices)
			write_debug_matrices();

		calculate_next_testvectors();

#ifdef UG_PARALLEL
		IF_DEBUG(LIB_ALG_AMG, 3)
		{
			UG_LOG("AH.get_master_layout:\n");
			PrintLayout(AH.get_master_layout());

			UG_LOG("AH.get_slave_layout :\n");
			PrintLayout(AH.get_slave_layout());
			UG_LOG("AH Layouts :\n");
			PrintLayout(A_OL2.get_communicator(), AH.get_master_layout(), AH.get_slave_layout());
		}
#endif

	}

private:
	//	get_aggressive_coarsening_interpolation
	//! tries for all nodes which are uninterpolateable an indirect interpolation
	//! \sa set_uninterpolateable_as_coarse
	void get_aggressive_coarsening_interpolation()
	{
		AMG_PROFILE_FUNC();

		for(size_t i=0; i<A.num_rows(); i++)
		{
			//UG_ASSERT(rating[i].is_valid_rating() == false, "node " << i << " has valid rating, but has not been processed yet?");
			if(rating[i].is_uninterpolateable() == false || rating.i_must_assign(i) == false) continue;

			// we have problems with nodes on the interface:
			// suppose i is a master node
			// 1. we would need the prolongation of the neighbors of i, which are slave0 nodes
			// 2. we could then indirectly interpolate from slave1 and master2 nodes
			// 3. on the slave side, we do not have the slave2 nodes, so we would have to add them somehow
			if(rating.is_inner_node(i))
				calculator.get_all_neighbors_interpolation(i, PoldIndices, rating);
			else
				rating.set_coarse(i);
		}
	}


	//	set_uninterpolateable_as_coarse
	//! sets all nodes which are uninterpolateable as coarse
	//! \sa get_aggressive_coarsening_interpolation
	void set_uninterpolateable_as_coarse()
	{
		AMG_PROFILE_FUNC();

		for(size_t i=0; i<rating.size(); i++)
		{
			if(rating[i].is_uninterpolateable() == false || rating.i_must_assign(i) == false) continue;
			rating.set_coarse(i);
		}
	}

	// create the bFine-Array for this level for F-smoothing
	void create_fine_marks()
	{
		AMG_PROFILE_FUNC();
		stopwatch SW;
		UG_DLOG(LIB_ALG_AMG, 1, std::endl << "create fine marks... "); if(bTiming) SW.start();
		stdvector<bool> &vFine = m_famg.levels[level]->is_fine;
		vFine.resize(A.num_rows());

		for(size_t i=0; i < A.num_rows(); i++)
			vFine[i] = rating[i].is_fine();
		if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms.");
	}

	//! - every coarse node c gets a unique new index newIndex[c]
	//! - calculate PnewIndices(r, newIndex[c]) = PoldIndex(r, c)
	//! - parentIndex[level+1][newIndex[i]] = i
	//! - next level master/slave-interfaces are updated with new indices.
	void create_new_indices()
	{
		AMG_PROFILE_BEGIN(create_new_indices);
		stopwatch SW;
		if(bTiming) { SW.start();UG_DLOG(LIB_ALG_AMG, 1, std::endl << "create new index... "); }

		size_t nrOfCoarse=0;
		stdvector<int> newIndex(overlapSize[1], -1);

		for(size_t r=0; r<A.num_rows(); r++)
		{
			for(typename prolongation_matrix_type::row_iterator it = PoldIndices.begin_row(r); it != PoldIndices.end_row(r); ++it)
			{
				size_t c = it.index();
				UG_ASSERT(rating[c].is_coarse(), c);
				if(newIndex[c] == -1)
					newIndex[c] = nrOfCoarse++;
			}
		}
		PnewIndices.resize(A.num_rows(), nrOfCoarse);
		for(size_t r=0; r<A.num_rows(); r++)
		{
			double maxCon=0;
			for(typename prolongation_matrix_type::row_iterator conn = PoldIndices.begin_row(r);
				conn != PoldIndices.end_row(r); ++conn)
				maxCon = BlockNorm(conn.value());


			for(typename prolongation_matrix_type::row_iterator conn = PoldIndices.begin_row(r);
					conn != PoldIndices.end_row(r); ++conn)
			{
				if(BlockNorm(conn.value()) < maxCon*m_famg.m_dEpsilonTr)
					continue;
				size_t c = conn.index();
				UG_ASSERT(rating[c].is_coarse(), "node " << c << " is not even coarse.");
				UG_ASSERT(newIndex[c] != -1, rating.info(c));
				// assure we are only interpolating Master0 and Slave0 nodes from Master0, Slave0 or Slave1 nodes.
				UG_ASSERT(rating.is_inner_node(r) || rating.is_master(r) || rating.is_slave(r, 0), rating.info(r));
				UG_ASSERT(rating.is_inner_node(c) || rating.is_master(c) || rating.is_slave(c, 0) || rating.is_slave(c, 1), rating.info(c));
				PnewIndices(r, newIndex[c]) = conn.value();
			}
		}

#ifdef UG_DEBUG
		for(size_t r=A.num_rows(); r<A_OL2.num_rows(); r++)
			UG_ASSERT(rating[r].is_fine() == false || rating[r].is_uncalculated_fine(), rating.info(r));
#endif
		UG_DLOG(LIB_ALG_AMG, 1, "rating.get_nr_of_coarse() = " << rating.get_nr_of_coarse() << ", nrOfCoarse = " << nrOfCoarse << "\n");
		PnewIndices.defragment();
#ifdef UG_PARALLEL
		PnewIndices.set_storage_type(PST_CONSISTENT);
#endif

		/*IF_DEBUG(LIB_ALG_AMG, 5)
		{
			UG_LOG("\n\nP:\n");
			PnewIndices.p();
			UG_LOG("\n\n");
		}*/

		if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms.");

		/// ---- create parent index ----
		AMG_PROFILE_NEXT(create_parent_index)
#ifndef UG_DEBUG
		//if(m_famg.m_writeMatrices)
		// we need this for restriction of test vectors
#endif
		{
			if(bTiming) { SW.start(); UG_DLOG(LIB_ALG_AMG, 1, std::endl << "create parentIndex... "); }

			stdvector<stdvector<int> > &parentIndex = m_famg.m_parentIndex;
			parentIndex.resize(level+2);
			parentIndex[level+1].resize(nrOfCoarse);
			for(size_t i=0; i < nrOfCoarse; i++) parentIndex[level+1][i] = -1;
			for(size_t i=0; i < overlapSize[1]; i++)
				if(newIndex[i] != -1)
					parentIndex[level+1][ newIndex[i] ] = i;

			if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms.");
		}

#ifdef UG_PARALLEL
		/// ---- update interfaces with new indices ----
		if(bTiming) { SW.start();UG_DLOG(LIB_ALG_AMG, 1, std::endl << "update interfaces with new indices... "); }
		AMG_PROFILE_NEXT(update_interfaces_with_new_indices)
		ReplaceIndicesInLayout(AH.get_master_layout(), newIndex);
		ReplaceIndicesInLayout(AH.get_slave_layout(), newIndex);
		if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms.\n");
		//PrintLayout(A_OL2.get_communicator(), nextLevelMasterLayout, nextLevelSlaveLayout);
#endif

		// create interface with newIndex
	}

	void write_debug_matrices();
	void write_debug_matrix_markers();
	template<typename TMatrix>
	void write_debug_matrix(TMatrix &mat, size_t fromlevel, size_t tolevel, const char *name);


	void on_demand_coarsening();
	void precalculate_coarsening();
	void calculate_all_possible_parent_pairs();

	void calculate_testvectors();
	void calculate_next_testvectors();

#ifdef UG_PARALLEL
	void color_process_graph();
	void receive_coarsening_from_processes_with_lower_color();
	void send_coarsening_data_to_processes_with_higher_color();
	void create_OL2_matrix();
	void add_connections_between_slave_nodes(IndexLayout &masterLayout, IndexLayout slaveLayout);
	void communicate_prolongation(); // in famg_communicate_prolongation.h
	void calculate_uncalculated_fine_nodes();
	void update_interface_with_newIndex(IndexLayout &layout, IndexLayout &nextLevelLayout, stdvector<int> &newIndex);
	IndexLayout OL1MasterLayout, OL1SlaveLayout;
	IndexLayout OL2MasterLayout, OL2SlaveLayout;
#endif


	std::vector<size_t> overlapSize;
};


template<>
void FAMG<CPUAlgebra>::get_testvectors(stdvector<vector_type> &testvectors, stdvector<double> &omega)
{
	AMG_PROFILE_FUNC();
	testvectors.resize(m_vVectorWriters.size() + m_testvectors.size());
	omega.resize(m_vVectorWriters.size() + m_testvectors.size());

	for(size_t i=0; i<m_testvectors.size(); i++)
	{
		testvectors[i] = m_testvectors[i];
		omega[i] = m_omegaVectors[i];
	}

	for(size_t i=0; i<m_vVectorWriters.size(); i++)
	{
		size_t index = i+m_testvectors.size();
		m_vVectorWriters[i]->update(testvectors[index]);
		omega[index] = m_omegaVectorWriters[i];

	}
}

#ifdef UG_PARALLEL
void eh( MPI_Comm *comm, int *err, ... )
{
	UG_LOG("MPI ERROR!\n");

	char error_string[256];
	int length_of_error_string=256, error_class;

	MPI_Error_class(*err, &error_class);
	MPI_Error_string(error_class, error_string, &length_of_error_string);
	UG_ASSERT(0,"MPI Error: " << error_string << "\n" );
	return;
}
#endif


template<>
void FAMG<CPUAlgebra>::c_create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
		prolongation_matrix_type &P, size_t level)
{
#ifdef UG_PARALLEL
    MPI_Errhandler newerr;

	MPI_Comm_create_errhandler( eh, &newerr );
	MPI_Comm_set_errhandler( MPI_COMM_WORLD, newerr );
	  //MPI_Comm_call_errhandler( MPI_COMM_WORLD, MPI_ERR_OTHER );
#endif


	if(level == 0)
	{
		testvectors.clear();
		get_testvectors(testvectors, omega);
	}

	if(m_writeMatrices && m_writeTestvectors)
	{
		for(size_t i=0; i<testvectors.size(); i++)
			WriteVectorToConnectionViewer(GetProcFilename(m_writeMatrixPath, ToString("testvector") + ToString(i) + ToString("_L") + ToString(level), ".vec").c_str(),
					testvectors[i], &m_amghelper.positions[level][0], m_dbgDimension);
	}

	AMG_PROFILE_FUNC();
	//UG_ASSERT(testvectors.size() > 0, "we need at least one testvector.");

	// testvectors will be altered by FAMGLevelCalculator

	// UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, 4);
	FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type> dummy(*this, AH, R, A, P,
			level, testvectors, omega);
	dummy.do_calculation();

}


} // namespace ug

#ifdef UG_PARALLEL
#include "famg_parallel_coarsening_impl.h"
#include "famg_communicate_prolongation.h"
#endif

#include "famg_debug_impl.h"
#include "famg_on_demand_coarsening_impl.h"
#include "famg_precalculate_coarsening_impl.h"
#include "famg_testvectors.h"
#endif //  __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__
