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

#include "famg.h"
#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/common/stl_debug.h"
#include "amg_debug_helper.h"

#include "stopwatch.h"
#include "common/assert.h"
#include "maxheap.h"
#include "famg.h"
#include <set>
#include <string>
#include <stack>

std::stack<int> g_DebugLevelStack;

#define PUSH_DEBUG_LEVEL() g_DebugLevelStack.push(ug::GetLogAssistant().get_debug_level(LIB_ALG_AMG))
#define SET_AND_PUSH_DEBUG_LEVEL(level) g_DebugLevelStack.push(ug::GetLogAssistant().get_debug_level(LIB_ALG_AMG)); ug::GetLogAssistant().set_debug_level(level)
#define POP_DEBUG_LEVEL() ug::GetLogAssistant().set_debug_level(g_DebugLevelStack.pop())


#define SET_DEBUG_LEVEL_OVERLAP() UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, 0); UG_SET_DEBUG_LEVEL(LIB_ALG_MATRIX, 0)
#define SET_DEBUG_LEVEL_TESTVECTOR_CALC() UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, 1); UG_SET_DEBUG_LEVEL(LIB_ALG_MATRIX, 0)
#define SET_DEBUG_LEVEL_PHASE_3() UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, 1)
#define SET_DEBUG_LEVEL_CALCULATE_PARENT_PAIRS() UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, 1)
#define SET_DEBUG_LEVEL_COLORING() UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, 1)
#define SET_DEBUG_LEVEL_RECV_COARSENING() UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, 1)
#define SET_DEBUG_LEVEL_GET_RATINGS() UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, 1)
#define SET_DEBUG_LEVEL_PRECALCULATE_COARSENING() UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, 1)
#define SET_DEBUG_LEVEL_AGGRESSIVE_COARSENING() UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, 1)
#define SET_DEBUG_LEVEL_SEND_COARSENING() UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, 1)
#define SET_DEBUG_LEVEL_COMMUNICATE_PROLONGATION() UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, 1)
#define SET_DEBUG_LEVEL_AFTER_COMMUNICATE_PROLONGATION() UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, 1)

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


#define AMG_WRITE_MATRICES_MAX (200*200)
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
#include "famg_testvectors.h"
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
	famg<CPUAlgebra> &m_famg;

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

public:
	// todo: clean up this mess of constructor
	FAMGLevelCalculator(famg<CPUAlgebra> &f,
			matrix_type &_AH, prolongation_matrix_type &_R,  const matrix_type &_A,
			prolongation_matrix_type &_P, size_t _level,
			stdvector< vector_type > &testvectors, stdvector<double> &omega)
	: m_famg(f), AH(_AH), A(_A), R(_R), PnewIndices(_P), level(_level),
			calculator(A, A_OL2, m_famg.get_delta(), m_famg.get_theta(),
					m_famg.get_damping_for_smoother_in_interpolation_calculation(), testvectors, omega),
			rating(PoldIndices, _level, m_famg.m_amghelper),
			m_testvectors(testvectors)
#ifndef UG_PARALLEL
			, A_OL2(A)
#else
			// TODO: fix that, its a memory leak
			, nextLevelMasterLayout(*(new IndexLayout)), nextLevelSlaveLayout(*(new IndexLayout))
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


		SET_DEBUG_LEVEL_OVERLAP();
		// 1. Overlap calculation
		//---------------------------
#ifdef UG_PARALLEL
		create_OL2_matrix();
#else
		rating.create(A_OL2.num_rows());
#endif
		size_t N = A_OL2.num_rows();

		SET_DEBUG_LEVEL_TESTVECTOR_CALC();
		// 2. global Testvector calculation (damping)
		//-----------------------------------------------
		// todo: all global?
		UG_DLOG(LIB_ALG_AMG, 1, "\ncalculating testvector... ");
		if(bTiming) SW.start();

		for(size_t i=0; i<m_testvectors.size(); i++)
		{
#ifdef UG_PARALLEL
			m_testvectors[i].set_storage_type(PST_CONSISTENT);
#endif
			CalculateTestvector(A_OL2,
					m_testvectors[i], m_famg.get_testvector_damps());
		}

		if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

		// 3. create heap, P, SymmNeighGraph
		//-------------------------------------
		SET_DEBUG_LEVEL_PHASE_3();
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

		bool bUsePrecalculate = true;
		if(bUsePrecalculate)
		{
			// get possible parent nodes
			UG_DLOG(LIB_ALG_AMG, 1, "\ncreating possible parent list... "); if(bTiming) SW.start();
			SET_DEBUG_LEVEL_CALCULATE_PARENT_PAIRS();
			calculate_all_possible_parent_pairs();

			if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");
		}


	#ifdef UG_PARALLEL
		// 4. coloring in parallel
		//-------------------------------------
		SET_DEBUG_LEVEL_COLORING();
		color_process_graph();
		SET_DEBUG_LEVEL_RECV_COARSENING();
		receive_coarsening_from_processes_with_lower_color();
	#endif

		rating.calculate_unassigned();

		// 5. do coarsening
		//-------------------------
		if(bUsePrecalculate)
		{
			// calculate ratings (not precalculateable because of coarse/uninterpolateable)
			SET_DEBUG_LEVEL_GET_RATINGS();
			UG_DLOG(LIB_ALG_AMG, 1, std::endl << "calculate ratings... "); if(bTiming) SW.start();
			GetRatings(possible_parents, rating, heap);
			if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

			//heap.print();

			// do coarsening
			SET_DEBUG_LEVEL_PRECALCULATE_COARSENING();
			UG_DLOG(LIB_ALG_AMG, 1, std::endl << "coarsening... "); if(bTiming) SW.start();
			precalculate_coarsening();
			if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms.");
		}
		else
		{
			UG_DLOG(LIB_ALG_AMG, 1, std::endl << "other coarsening... "); if(bTiming) SW.start();
			on_demand_coarsening();
			if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms.");
		}


		// 7. do aggressive coarsening
		//-----------------------------------------------
		// this has to be AFTER communication of coarse/fine markers,
		// since on border nodes, we will have

		SET_DEBUG_LEVEL_AGGRESSIVE_COARSENING();
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
			SET_DEBUG_LEVEL_SEND_COARSENING();
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
			SET_DEBUG_LEVEL_COMMUNICATE_PROLONGATION();
			communicate_prolongation();
		#endif

		// [ debug output
		write_debug_matrix_markers();
		// ]

		SET_DEBUG_LEVEL_AFTER_COMMUNICATE_PROLONGATION();
		IF_DEBUG(LIB_ALG_AMG, 3)
		{	rating.print(); 	}

		create_new_index();
		PnewIndices.defragment();

	#ifdef UG_PARALLEL
		PnewIndices.set_storage_type(PST_CONSISTENT);
	#endif

		UG_DLOG(LIB_ALG_AMG, 1, std::endl << "create fine marks... "); if(bTiming) SW.start();
		create_fine_marks();
		if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms.");
		//UG_DLOG(LIB_ALG_AMG, 1, "parentIndex level " << level << "\n")
		//for(size_t i=0; i<rating.get_nr_of_coarse(); i++) { UG_ASSERT(amghelper.parentIndex[level+1][i] != -1, i << " == -1???"); UG_LOG(i << " = " << amghelper.parentIndex[level+1][i] << "\n"); }


		IF_DEBUG(LIB_ALG_AMG, 4)
		{
			UG_LOG("\n\nP:\n");
			PnewIndices.p();
			UG_LOG("\n\n");
		}

		UG_DLOG(LIB_ALG_AMG, 1, std::endl << rating.get_nr_of_coarse() << " coarse, " << N - rating.get_unassigned() - rating.get_nr_of_coarse() << " fine.");

		// 8. construct restriction R = I_{h->2h}
		//-----------------------------------------

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

		// 9. create Galerkin product AH = R A P
		//-----------------------------------------

		UG_DLOG(LIB_ALG_AMG, 1, "\ngalerkin product... "); if(bTiming) SW.start();

		// AH = R A P
		PROFILE_BEGIN(FAMGCreateAsMultiplyOf)
			CreateAsMultiplyOf(AH, R, A, PnewIndices, 0.01);
		PROFILE_END();

		if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

		// finalize
		if(bTiming) { UG_DLOG(LIB_ALG_AMG, 1, std::endl << "Finalizing.."); SW.start(); }
		AH.defragment();
		if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, " took " << SW.ms() << " ms\n");

	#ifdef UG_PARALLEL
		AH.set_storage_type(PST_ADDITIVE);
		AH.set_layouts(nextLevelMasterLayout, nextLevelSlaveLayout);
	#endif


	#ifdef FAMG_PRINT_AH
		IF_DEBUG(LIB_ALG_AMG, 3)
		{
			UG_DLOG(LIB_ALG_AMG, 1, "AH level " << level << std::endl);
			AH.print();
		}
	#endif

		// 10.
		//-----------------------------------------

		write_debug_matrices();

		// todo: remove dynamic cast, change big_testvector to parallel
		for(size_t i=0; i<m_testvectors.size(); i++)
			CalculateNextTestvector(R, m_testvectors[i]);

#ifdef UG_PARALLEL
		IF_DEBUG(LIB_ALG_AMG, 3)
		{
			UG_LOG("nextLevelMasterLayout :\n");
			PrintLayout(nextLevelMasterLayout);

			UG_LOG("nextLevelSlaveLayout :\n");
			PrintLayout(nextLevelSlaveLayout);
			UG_LOG("nextLevelLayout :\n");
			PrintLayout(A_OL2.get_communicator(), nextLevelMasterLayout, nextLevelSlaveLayout);
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

			calculator.get_all_neighbors_interpolation(i, PoldIndices, rating);
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

	void create_fine_marks()
	{
		AMG_PROFILE_FUNC();
		m_famg.is_fine.resize(level+1);
		stdvector<bool> &vFine = m_famg.is_fine[level];
		vFine.resize(A.num_rows());

		for(size_t i=0; i < A.num_rows(); i++)
			vFine[i] = rating[i].is_fine();
	}

	//! creates the parentIndex struct for display
	void create_new_index()
	{
		AMG_PROFILE_BEGIN(create_new_index);
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

		for(size_t r=A.num_rows(); r<A_OL2.num_rows(); r++)
			UG_ASSERT(rating[r].is_fine() == false || rating[r].is_uncalculated_fine(), rating.info(r));

		PnewIndices.resize(A.num_rows(), nrOfCoarse);
		UG_DLOG(LIB_ALG_AMG, 1, "rating.get_nr_of_coarse() = " << rating.get_nr_of_coarse() << ", nrOfCoarse = " << nrOfCoarse << "\n");

		if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms.");

		AMG_PROFILE_NEXT(create_parent_index)
		if(m_famg.m_writeMatrices)
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
		if(bTiming) { SW.start();UG_DLOG(LIB_ALG_AMG, 1, std::endl << "create interfaces... "); }
		AMG_PROFILE_NEXT(create_interfaces)
		update_interface_with_newIndex(A.get_master_layout(), nextLevelMasterLayout, newIndex);
		update_interface_with_newIndex(A.get_slave_layout(), nextLevelSlaveLayout, newIndex);
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
void famg<CPUAlgebra>::get_testvectors(stdvector<vector_type> &testvectors, stdvector<double> &omega)
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


		if(m_writeMatrices && m_writeTestvectors)
		{
			vector_type &vec = testvectors[index];
			std::fstream file((m_writeMatrixPath + "testvector" + ToString(i) + ".values").c_str(),
					std::ios::out);
			for(size_t j=0; j<vec.size(); j++)
				file << j << " " << (vec[j]) << std::endl;


			WriteVectorToConnectionViewer((m_writeMatrixPath + "testvector" + ToString(i) +".mat").c_str(),
				vec, &m_amghelper.positions[0], 2);
		}

	}
}


template<>
void famg<CPUAlgebra>::c_create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
		prolongation_matrix_type &P, size_t level)
{
	AMG_PROFILE_FUNC();
	stdvector< vector_type > testvectors;
	stdvector<double> omega;
	get_testvectors(testvectors, omega);

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
#endif //  __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__
