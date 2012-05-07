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
		stopwatch SW;
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
		stopwatch SW;
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
		stopwatch SW, SWwhole; SWwhole.start();

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
			WriteVectorToConnectionViewer(
				m_famg.m_writeMatrixPath + std::string("AMG_fvalues_L") + ToString(level) + ".vec",
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

template<>
void FAMG<CPUAlgebra>::get_testvectors_from_matrix_rows
	(const matrix_type &A, stdvector<vector_type> &testvectors, stdvector<double> &omega)
{
	AMG_PROFILE_FUNC();

	testvectors.resize(1);
	omega.resize(1);
	vector_type &v = testvectors[0];
	v.resize(A.num_rows());
	for(size_t i=0; i<A.num_rows(); i++)
	{
		if(A.is_isolated(i))
			v[i] = 0.0;
		else
			v[i] = 1.0;
	}
}


template<>
void FAMG<CPUAlgebra>::precalc_level(size_t level)
{

	AMG_PROFILE_FUNC();
	//UG_LOG("\n\n\nprecalc!\n\n\n");
	UG_ASSERT(m_testvectorsmoother != NULL, "please provide a testvector smoother.");


	if(level == 0)
	{
		testvectors.clear();
		if(m_bTestvectorsFromMatrixRows)
			get_testvectors_from_matrix_rows(*levels[level]->pAgglomeratedA, testvectors, omega);
		else
			get_testvectors(testvectors, omega);
	}

#ifdef UG_PARALLEL
	if((level != 0 || m_bTestvectorsFromMatrixRows == false)
			&& levels[level]->bHasBeenMerged)
	{
		vector_type t;
		if(AMGBase<CPUAlgebra>::isMergingMaster(level))
			t.resize(levels[level]->pAgglomeratedA->num_rows());
		for(size_t i=0; i<testvectors.size(); i++)
		{
			t.set_storage_type(PST_CONSISTENT);
			testvectors[i].set_storage_type(PST_CONSISTENT);
			gather_vertical(testvectors[i], t, level, PST_CONSISTENT);
			testvectors[i] = t;
		}
	}
#endif
}

template<>
void FAMG<CPUAlgebra>::c_create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
		prolongation_matrix_type &P, size_t level)
{
	AMG_PROFILE_FUNC();

	UG_ASSERT(m_testvectorsmoother != NULL, "please provide a testvector smoother.");

	UG_ASSERT(testvectors.size() != 0, "testvectors?");
	UG_ASSERT(testvectors[0].size() == A.num_rows(), "testvectorsize = " << testvectors[0].size() << " != A.num_rows() = " << A.num_rows() << " ?");

	if(m_writeMatrices && m_writeTestvectors)
	{
		for(size_t i=0; i<testvectors.size(); i++)
			WriteVectorToConnectionViewer((m_writeMatrixPath + ToString("testvector") + ToString(i) + ToString("_L") + ToString(level) + ".vec").c_str(),
					testvectors[i], &m_amghelper.positions[level][0], m_dbgDimension);
	}

	//UG_ASSERT(testvectors.size() > 0, "we need at least one testvector.");

	// testvectors will be altered by FAMGLevelCalculator

	// UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, 4);

	if(m_dPrereduceAToStrongParameter > 0.0)
	{
		matrix_type Aeps;
		ReduceToStrongConnections(Aeps, A, m_dPrereduceAToStrongParameter);

		FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type> dummy(*this, AH, R, Aeps, A, P, level, testvectors, omega);
		dummy.do_calculation();
	}
	else
	{
		FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type> dummy(*this, AH, R, A, A, P, level, testvectors, omega);
		dummy.do_calculation();
	}

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//       check_testvector
//--------------------------------
/**
 */
template<>
bool FAMG<CPUAlgebra>::check_testvector()
{
	AMG_PROFILE_FUNC();

	UG_LOG("\n");
	UG_LOG("            check_testvector\n");
	UG_LOG("==========================================\n")
	UG_LOG("check if testvectors are interpolated exactly\n")
	if (0){
		int level = AMGBase<CPUAlgebra>::m_usedLevels-1;
		matrix_type &A = *AMGBase<CPUAlgebra>::levels[level]->pA;
		for(size_t i=0; i<A.num_rows(); i++)
		{
			double sum=0;
			for(matrix_type::row_iterator it = A.begin_row(i); it != A.end_row(i); ++it)
				sum += it.value();
			if(dabs(sum) > 0.1)
			{
				UG_LOG("Row " << i << " has sum " << sum << ": ");
				A.pr(i);
				UG_LOG("\n");
			}
		}
	}

	size_t preSmooth = AMGBase<CPUAlgebra>::get_num_presmooth();
	size_t postSmooth = AMGBase<CPUAlgebra>::get_num_postsmooth();
	AMGBase<CPUAlgebra>::set_num_presmooth(get_testvector_smooths());
	AMGBase<CPUAlgebra>::set_num_postsmooth(get_testvector_smooths());






	for(size_t level=0; ; level++)
	{
		// get testvectors, get agglomerated tv
		precalc_level(level);
#ifdef UG_PARALLEL
		if(AMGBase<CPUAlgebra>::isMergingSlave(level))
		{
			UG_LOG("merged.\n");
			break;
		}
#endif
		UG_LOG("Check Level " << level << "\n-----------------------------------\n")

		SmartPtr<matrix_type> pA, pAH;
		pA = AMGBase<CPUAlgebra>::levels[level]->pAgglomeratedA;
		pAH = AMGBase<CPUAlgebra>::levels[level+1]->pAgglomeratedA;

		matrix_type &A = *pA;
		matrix_type &AH = *pAH;
		matrix_type &R = AMGBase<CPUAlgebra>::levels[level]->R;
		matrix_type &P = AMGBase<CPUAlgebra>::levels[level]->P;
		vector_type c, d, tv;

	#ifdef UG_PARALLEL
		SetParallelVectorAsMatrix(c, A, PST_CONSISTENT);
		SetParallelVectorAsMatrix(d, A, PST_CONSISTENT);
		SetParallelVectorAsMatrix(tv, A, PST_CONSISTENT);
	#endif
		tv.resize(A.num_rows());
		d.resize(A.num_rows());
		c.resize(A.num_rows());
		for(size_t i=0; i<d.size(); i++)
			tv[i] = testvectors[0][i];
		A.apply(d, tv);
		c.set(0.0);
		checkResult res;

		/*PRINTLAYOUT(A.get_process_communicator(), A.get_communicator(), d.get_master_layout(), d.get_slave_layout());
		d.change_storage_type(PST_CONSISTENT);

		UG_LOG("\n\ntv level " << level << "\n");
		tv.print();
		UG_LOG("\n\nA level " << level << "\n");
		A.print();
		UG_LOG("\n\ntestvectors level " << level << "\n");
		testvectors[0].print();
		UG_LOG("\n\nd level " << level << "\n");
		d.print();

		d.change_storage_type(PST_ADDITIVE);*/

		AMGBase<CPUAlgebra>::check_level(c, d, A, level, res, &tv);

		if(level+1 < AMGBase<CPUAlgebra>::m_usedLevels-1)
		{
			matrix_type Aeps;
			ReduceToStrongConnections(Aeps, A, m_dPrereduceAToStrongParameter);
			FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type> dummy(*this, AH, R, Aeps, A, P, level, testvectors, omega);
			dummy.onlyTV();
		}
		else break;

	}
	AMGBase<CPUAlgebra>::set_num_presmooth(preSmooth);
	AMGBase<CPUAlgebra>::set_num_postsmooth(postSmooth);
	return true;
}
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
