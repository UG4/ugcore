/**
 * \file famg_level_calculator.h
 *
 * \author Martin Rupp
 *
 * \date 7.12.2010
 *
 * implementation file for famg
 *
 * Goethe-Center for Scientific Computing 2010-2012.
 *
 */


#ifndef __H__LIB_ALGEBRA__AMG__FAMG_LEVEL_CALCULATOR_H__
#define __H__LIB_ALGEBRA__AMG__FAMG_LEVEL_CALCULATORH__

//#include "sparsematrix_util.h"

#include "common/assert.h"
#include "common/log.h"

#include "famg.h"
#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/cpu_algebra_types.h"
#include "lib_algebra/common/stl_debug.h"
#include "../amg_debug_helper.h"

#include "common/stopwatch.h"
#include "common/assert.h"
#include "lib_algebra/common/heap/maxheap.h"
#include "famg.h"
#include <set>
#include <string>
#include <stack>
#include "../amg_misc.h"

#ifdef UG_PARALLEL
#include "lib_algebra/parallelization/parallel_matrix_overlap_impl.h"
#include "lib_algebra/parallelization/communication_scheme.h"
#include "pcl/pcl_layout_tests.h"
#endif

#include "lib_algebra/common/connection_viewer_output.h"

//#include "../row_sender.h"



std::stack<int> g_DebugLevelStack;


#ifdef UG_DEBUG
#define GET_DEBUG_LEVEL(tag) ug::GetLogAssistant().get_debug_level(ug::LogAssistant::tag)
#else
#define GET_DEBUG_LEVEL(tag) 0
#endif
#define PUSH_DEBUG_LEVEL() g_DebugLevelStack.push(GET_DEBUG_LEVEL(LIB_ALG_AMG))
#define SET_AND_PUSH_DEBUG_LEVEL(level) g_DebugLevelStack.push(GET_DEBUG_LEVEL(LIB_ALG_AMG)); ug::GetLogAssistant().set_debug_level(level)
#define POP_DEBUG_LEVEL() ug::GetLogAssistant().set_debug_level(g_DebugLevelStack.pop())

/*
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
}*/


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

namespace ug{

#ifdef UG_PARALLEL
	void MyPrintLayout(pcl::ProcessCommunicator &pc, pcl::InterfaceCommunicator<IndexLayout> &communicator, IndexLayout &layout1, IndexLayout &layout2, const char *name1, const char *name2)
{
	UG_LOG("\n========================================\n");
	UG_LOG(name1<<"\n");
	PrintLayout(layout1);
	UG_LOG("\n"<<name2<<"\n");
	PrintLayout(layout2);
	UG_LOG("\n"<<name1 << "-" << name2 <<"\n");
	PrintLayout(pc, communicator, layout1, layout2);
	UG_LOG("\n\n");
}
#endif




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
static void CreateSymmConnectivityGraph(const matrix_type &A, cgraph &SymmNeighGraph)
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

	// the fine grid matrix (const!)
	const matrix_type &Aorig;

	// interpolation/restriction matrices to be filled
	prolongation_matrix_type &R, &PnewIndices;

	size_t level; // current amg level

#ifdef UG_PARALLEL
	ParallelNodes PN;
	//! layout for sending/receiving coarsening information
	//! note: this is not a classical master/slave distribution!
	IndexLayout OLCoarseningReceiveLayout, OLCoarseningSendLayout;

	//! processes with higher/lower color to distribute coarsening information
	std::vector<int> processesWithLowerColor, processesWithHigherColor;
	int m_myColor;

	// overlap 2 matrix in parallel case
	matrix_type A_OL2;
	IndexLayout &nextLevelMasterLayout;
	IndexLayout &nextLevelSlaveLayout;
	IndexLayout OL1MasterLayout, OL1SlaveLayout;
	IndexLayout OL2MasterLayout, OL2SlaveLayout;
	std::vector<size_t> overlapSize;
	std::vector<IndexLayout> masterLayouts;
	std::vector<IndexLayout> slaveLayouts;
#else
	const matrix_type &A_OL2;
#endif
	FAMGNodes rating;
	stdvector< vector_type > &m_testvectors;

	// the calculator which calculates our interpolation weighs
	FAMGInterpolationCalculator<matrix_type, vector_type> calculator;

	bool bTiming; //< if true, timing will be outputed

	// this structure holds coarse/fine information as well as
	// master/slave. it tries to replace the missing "node" object here


	//! list of possible interpolating parents for each node
	stdvector<stdvector<neighborstruct2> > possible_parents;

	stdvector<bool> prolongation_calculated;

	//! heap used for sorting of ratings
	maxheap<FAMGNode> heap;

	//! used to determine which neighbors' rating needs to be updated:
	cgraph SymmNeighGraph;

	stdvector<double> fvalues;
	// our testvectors


	prolongation_matrix_type PoldIndices;

	void create_restriction()
	{
		AMG_PROFILE_FUNC();
		Stopwatch SW;
		UG_DLOG(LIB_ALG_AMG, 0, std::endl << "restriction... "); if(bTiming) SW.start();
		// construct restriction R = I_{h -> 2h}
		R.set_as_transpose_of(PnewIndices);
		// R is already defragmented
		if(bTiming) UG_DLOG(LIB_ALG_AMG, 0, "took " << SW.ms() << " ms.");

	#ifdef UG_PARALLEL
		R.set_storage_type(PST_CONSISTENT);
	#endif

	#ifdef FAMG_PRINT_R
		IF_DEBUG(LIB_ALG_AMG, 4)
		{
			UG_DLOG(LIB_ALG_AMG, 0, std::endl << "Restriction level " << level << std::endl);
			R.print();
		}
	#endif
	}

	void create_galerkin_product()
	{
		AMG_PROFILE_FUNC();
		Stopwatch SW;
		UG_DLOG(LIB_ALG_AMG, 0, "\ngalerkin product... "); if(bTiming) SW.start();
		// AH = R A P

		/* matrix_type AH2;
		matrix_type RoldIndices;
		RoldIndices.set_as_transpose_of(PoldIndices);
		CreateAsMultiplyOf2(AH2, RoldIndices, A, PoldIndices); */

		// R.print("R");
		// A.print("A");
		// PnewIndices.print("P");
		CreateAsMultiplyOf(AH, R, Aorig, PnewIndices, m_famg.get_galerkin_truncation());
		// AH.print();
		if(bTiming) UG_DLOG(LIB_ALG_AMG, 0, "took " << SW.ms() << " ms");

		// finalize
		if(bTiming) { UG_DLOG(LIB_ALG_AMG, 0, std::endl << "Finalizing.."); SW.start(); }
		AH.defragment();
		if(bTiming) UG_DLOG(LIB_ALG_AMG, 0, " took " << SW.ms() << " ms\n");



	#ifdef UG_PARALLEL
		AH.set_storage_type(PST_ADDITIVE);
		AH.set_communicator(A_OL2.communicator());
		AH.set_process_communicator(A_OL2.process_communicator());
		// AH.set_layouts(nextLevelMasterLayout, nextLevelSlaveLayout); // is already set
	#endif


		PROFILE_END();

	#ifdef FAMG_PRINT_AH
		IF_DEBUG(LIB_ALG_AMG, 3)
		{
			UG_DLOG(LIB_ALG_AMG, 0, "AH level " << level << std::endl);
			AH.print();
		}

	#endif
	}
public:
	// todo: clean up this mess of constructor
	FAMGLevelCalculator(FAMG<CPUAlgebra> &f,
			matrix_type &_AH, prolongation_matrix_type &_R,  const matrix_type &_A, const matrix_type &_Aorig,
			prolongation_matrix_type &_P, size_t _level,
			stdvector< vector_type > &testvectors, stdvector<double> &omega)
	: m_famg(f), AH(_AH), A(_A), Aorig(_Aorig), R(_R), PnewIndices(_P), level(_level),
#ifdef UG_PARALLEL

			PN(const_cast<matrix_type&>(A).communicator(), const_cast<matrix_type&>(A).master_layout(),
						const_cast<matrix_type&>(A).slave_layout(), A.num_rows()),
			nextLevelMasterLayout(AH.master_layout()), nextLevelSlaveLayout(AH.slave_layout()),
			rating(PoldIndices, _level, m_famg.m_amghelper, PN),
#else
			A_OL2(A),
			rating(PoldIndices, _level, m_famg.m_amghelper),
#endif
			m_testvectors(testvectors),
			calculator(A, A_OL2, m_famg, testvectors, omega, fvalues)
	{
	}

public:
	void do_calculation()
	{
		AMG_PROFILE_FUNC();

		UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, 1);
		UG_DLOG(LIB_ALG_AMG, 0, "Creating level " << level << ". (" << A.num_rows() << " nodes)" << std::fixed);
		bTiming=true;
		Stopwatch SW, SWwhole; SWwhole.start();

		// 1. Overlap calculation
		//---------------------------
#ifdef UG_PARALLEL
		create_OL2_matrix();
#else
		rating.create(A_OL2.num_rows());
#endif
		size_t N = A_OL2.num_rows();

		if(m_famg.m_bWriteFValues)	fvalues.resize(N, 0.0);
		else 						fvalues.clear();

		// 2. global Testvector calculation (damping)
		//-----------------------------------------------

		calculate_testvectors();
		get_testvectors_on_OL2();
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

			UG_DLOG(LIB_ALG_AMG, 0, "\nCreate P, SymmNeighGraph... "); if(bTiming) SW.start();
			PoldIndices.resize(N, N);
			// get neighboring information
			SymmNeighGraph.resize(N);
			CreateSymmConnectivityGraph(A_OL2, SymmNeighGraph);
			if(bTiming) UG_DLOG(LIB_ALG_AMG, 0, "took " << SW.ms() << " ms");

			if(bUsePrecalculate)
				// get possible parent nodes
				calculate_all_possible_parent_pairs();

#ifdef UG_PARALLEL
			// 4. coloring in parallel
			//-------------------------------------
			color_process_graph();
			receive_coarsening_from_processes_with_lower_color();
#endif
			rating.calculate_unassigned();

			// 5. do coarsening
			//-------------------------
			if(bUsePrecalculate)
			{
				// calculate ratings (not precalculateable because of coarse/uninterpolateable)
				UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, m_famg.iDebugLevelGetRatings);
				UG_DLOG(LIB_ALG_AMG, 0, std::endl << "calculate ratings... "); if(bTiming) SW.start();
				GetRatings(possible_parents, rating, heap);
				if(bTiming) UG_DLOG(LIB_ALG_AMG, 0, "took " << SW.ms() << " ms");

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
			UG_DLOG(LIB_ALG_AMG, 0, std::endl << "second coarsening... "); if(bTiming) SW.start();
			if(m_famg.get_aggressive_coarsening() == true)
				get_aggressive_coarsening_interpolation();
			else
				set_uninterpolateable_as_coarse();
			if(bTiming) UG_DLOG(LIB_ALG_AMG, 0, "took " << SW.ms() << " ms.");



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

			UG_DLOG(LIB_ALG_AMG, 0, std::endl << N - rating.get_unassigned() << " nodes assigned, " << rating.get_nr_of_coarse() << " coarse, "
					<< N - rating.get_unassigned() - rating.get_nr_of_coarse() << " fine, " << rating.get_unassigned() << " unassigned.");

			// set P to real size

			#ifdef UG_PARALLEL
			/*{
				pcl::InterfaceCommunicator<IndexLayout> &communicator = (const_cast<matrix_type&>(A)).communicator();

				ParallelNodes PN(OL2MasterLayout, OL2SlaveLayout, A_OL2.num_rows());

				communicate_prolongation(communicator, A.master_layout(), A.slave_layout(),
						nextLevelMasterLayout, nextLevelSlaveLayout, PN, rating);
			}*/
			#endif
		}
		else
		{
			if(m_famg.m_writeMatrices)	m_famg.write_debug_matrix_markers(level, rating);
			PoldIndices.resize(N, N);
			rs_amg_external_coarsening();
			external_coarsening_calculate_prolongation();
		}



		UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, m_famg.iDebugLevelAfterCommunicateProlongation);
		IF_DEBUG(LIB_ALG_AMG, 3)
		{	rating.print(); 	}


		create_new_indices();
		//for(size_t i=0; i<rating.get_nr_of_coarse(); i++) { UG_ASSERT(amghelper.parentIndex[level+1][i] != -1, i << " == -1???"); UG_LOG(i << " = " << amghelper.parentIndex[level+1][i] << "\n"); }

		//UG_DLOG(LIB_ALG_AMG, 0, std::endl << rating.get_nr_of_coarse() << " coarse.");

		// 8. construct restriction R = I_{h->2h}
		//-----------------------------------------
		create_restriction();


		for(size_t i=0; i<A.num_rows(); i++)
		{
			if(PoldIndices.num_connections(i) > 1 && (rating[i].is_fine() || rating[i].is_fine_direct()) == false)
			{	UG_LOG("level " << level << " node " << i << "\n");		}
		}
		// 9. create Galerkin product AH = R A P
		//-----------------------------------------

		create_galerkin_product();
		// 10.
		//-----------------------------------------
		// [ debug output
		if(m_famg.m_writeMatrices)
		{
			m_famg.write_debug_matrix_markers(level, rating);
			m_famg.write_debug_matrices(AH, R, A, PnewIndices, level);
		}
		// ]

		calculate_next_testvectors();

#ifdef UG_PARALLEL
		IF_DEBUG(LIB_ALG_AMG, 3)
		{
			PRINTLAYOUT(A_OL2.process_communicator(), A_OL2.communicator(), AH.master_layout(), AH.slave_layout());
		}
#endif

		if(m_famg.m_bWriteFValues)
		{
			std::string path=std::string("/level") + ToString(level) + "/";
			mkdir((std::string(m_famg.m_writeMatrixPath) + path).c_str(), 0777);
			WriteVectorToConnectionViewer(
				m_famg.m_writeMatrixPath + path + std::string("AMG_fvalues_L") + ToString(level) + ".vec",
				fvalues, &m_famg.m_amghelper.positions[level][0], 2);
		}
	}

private:

	void create_new_indices()
	{
		AMG_PROFILE_FUNC();
		UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, m_famg.iDebugLevelCommunicateProlongation);
		UG_SET_DEBUG_LEVEL(LIB_ALG_MATRIX, m_famg.iDebugLevelCommunicateProlongation);

#ifdef UG_PARALLEL
		PoldIndices.set_process_communicator(A.process_communicator());
		m_famg.parallel_process_prolongation(PoldIndices, PnewIndices, m_famg.m_dProlongationTruncation, level, rating,
				PN, false, nextLevelMasterLayout, nextLevelSlaveLayout);
		TESTLAYOUT(A_OL2.process_communicator(), A_OL2.communicator(), nextLevelMasterLayout, nextLevelSlaveLayout);
#else
		m_famg.serial_process_prolongation(PoldIndices, PnewIndices, m_famg.m_dProlongationTruncation, level, rating);
#endif

	}

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

	void on_demand_coarsening();
	void precalculate_coarsening();
	void calculate_all_possible_parent_pairs();

#ifdef UG_PARALLEL
	void color_process_graph();
	void receive_coarsening_from_processes_with_lower_color();
	void send_coarsening_data_to_processes_with_higher_color();
	void create_OL2_matrix();
	void add_connections_between_slave_nodes(IndexLayout &masterLayout, IndexLayout slaveLayout);
	void calculate_uncalculated_fine_nodes();
	void update_interface_with_newIndex(IndexLayout &layout, IndexLayout &nextLevelLayout, stdvector<int> &newIndex);
#endif
	void external_coarsening_calculate_prolongation();
	void rs_amg_external_coarsening();

	void calculate_testvectors();
	void get_testvectors_on_OL2();
	void calculate_next_testvectors();
public:
	void onlyTV()
	{
#ifdef UG_PARALLEL
		create_OL2_matrix();
#endif
		calculate_testvectors();
		get_testvectors_on_OL2();
		calculate_next_testvectors();
	}
};



} // namespace ug

#ifdef UG_PARALLEL
#include "famg_parallel_coarsening_impl.h"
#endif

#include "famg_debug_impl.h"
#include "famg_on_demand_coarsening_impl.h"
#include "famg_precalculate_coarsening_impl.h"
#include "famg_testvectors.h"
#include "famg_external_coarsening.h"
#endif //  __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__
