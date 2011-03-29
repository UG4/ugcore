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

//#define PROFILE_FAMG
#ifdef PROFILE_FAMG
	#define FAMG_PROFILE_FUNC()			PROFILE_FUNC()
	#define FAMG_PROFILE_BEGIN(name)	PROFILE_BEGIN(name)
	#define FAMG_PROFILE_END()			PROFILE_END()
#else
	#define FAMG_PROFILE_FUNC()
	#define FAMG_PROFILE_BEGIN(name)
	#define FAMG_PROFILE_END()
#endif

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


#ifdef UG_PARALLEL
std::string GetProcFilename(std::string name, std::string extension)
{
	return name + "_" + ToString(pcl::GetProcRank()) + extension;
}

#else
std::string GetProcFilename(string name, string extension)
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
	FAMG_PROFILE_FUNC();

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
	prolongation_matrix_type &R, &P;

	size_t level; // current amg level

	// the calculator which calculates our interpolation weighs
	FAMGInterpolationCalculator<matrix_type, vector_type> calculator;

	// this structure holds coarse/fine information as well as
	// master/slave. it tries to replace the missing "node" object here
	famg_nodes rating;

	//! list of possible interpolating parents for each node
	stdvector<stdvector<neighborstruct2> > possible_parents;

	//! heap used for sorting of ratings
	maxheap<famg_nodeinfo> heap;

	//! used to determine which neighbors' rating needs to be updated:
	cgraph SymmNeighGraph;

	// our testvectors
	stdvector< vector_type > m_testvectors;

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
	matrix_type &A_OL2;
#endif



public:
	// todo: clean up this mess of constructor
	FAMGLevelCalculator(famg<CPUAlgebra> &f,
			matrix_type &_AH, prolongation_matrix_type &_R,  const matrix_type &_A,
			prolongation_matrix_type &_P, size_t _level,
			stdvector< vector_type > &testvectors, stdvector<double> &omega)
	: m_famg(f), AH(_AH), A(_A), R(_R), P(_P), level(_level),
			calculator(A, A_OL2, m_famg.get_delta(), m_famg.get_theta(),
					m_famg.get_damping_for_smoother_in_interpolation_calculation(), testvectors, omega),
			rating(P, _level, m_famg.m_amghelper),
			m_testvectors(testvectors)
#ifndef UG_PARALLEL
			, A_OL2(A)
#else
			// TODO: fix that, its a memory leak
			, nextLevelMasterLayout(*(new IndexLayout)), nextLevelSlaveLayout(*(new IndexLayout))
#endif
	{
	}

private:
	//	get_aggressive_coarsening_interpolation
	//! tries for all nodes which are uninterpolateable an indirect interpolation
	//! \sa set_uninterpolateable_as_coarse
	void get_aggressive_coarsening_interpolation()
	{
		FAMG_PROFILE_FUNC();

		UG_SET_DEBUG_LEVELS(4);
		for(size_t i=0; i<A.num_rows(); i++)
		{
			UG_ASSERT(rating[i].is_valid_rating() == false, "node " << i << " has valid rating, but has not been processed yet?");
			if(rating[i].is_uninterpolateable() == false) continue;

			calculator.get_all_neighbors_interpolation(i, P, rating);
		}
	}


	//	set_uninterpolateable_as_coarse
	//! sets all nodes which are uninterpolateable as coarse
	//! \sa get_aggressive_coarsening_interpolation
	void set_uninterpolateable_as_coarse()
	{
		FAMG_PROFILE_FUNC();

		for(size_t i=0; i<rating.size(); i++)
		{
			UG_ASSERT(rating[i].is_valid_rating() == false, "node " << i << " has valid rating (neither coarse nor fine but interpolateable), but is not in heap anymore???");
			if(rating[i].is_uninterpolateable() == false) continue;
			rating.set_coarse(i);
		}
	}

	//! creates the parentIndex struct for display
	void create_parentIndex()
	{
		FAMG_PROFILE_FUNC();

		stdvector<stdvector<int> > &parentIndex = m_famg.m_parentIndex;
		parentIndex.resize(level+2);
		parentIndex[level+1].resize(rating.get_nr_of_coarse());
		for(size_t i=0; i < rating.get_nr_of_coarse(); i++) parentIndex[level+1][i] = -1;
		for(size_t i=0; i < rating.size(); i++) {
			if(rating.newIndex[i] != -1) parentIndex[level+1][ rating.newIndex[i] ] = i;
		}
	}


#ifdef UG_PARALLEL
	void color_process_graph();
	void receive_coarsening_from_processes_with_lower_color();
	void send_coarsening_data_to_processes_with_higher_color();
	void create_OL2_matrix();
	void add_connections_between_slave_nodes(IndexLayout &masterLayout, IndexLayout slaveLayout);
#endif

public:
	void do_calculation()
	{
		FAMG_PROFILE_FUNC();

		UG_LOG("Creating level " << level << ". (" << A.num_rows() << " nodes)" << std::fixed);
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

		// UG_SET_DEBUG_LEVELS(4);

		// 2. global Testvector calculation (damping)
		//-----------------------------------------------
		// todo: all global?
		UG_LOG("\ncalculating testvector... ");
		if(bTiming) SW.start();

		for(size_t i=0; i<m_testvectors.size(); i++)
		{
#ifdef UG_PARALLEL
			m_testvectors[i].set_storage_type(PST_CONSISTENT);
#endif
			CalculateTestvector(A_OL2,
					m_testvectors[i], m_famg.get_testvector_damps());
		}

		if(bTiming) UG_LOG("took " << SW.ms() << " ms");

		// 3. create heap, P, SymmNeighGraph
		//-------------------------------------
		heap.create(rating.nodes);

		IF_DEBUG(LIB_ALG_AMG, 3)
		{
			UG_LOG("\bThe Matrix A:\n\n");
			A.print();
			UG_LOG("\n\n");
		}

		UG_LOG("\nCreate P, SymmNeighGraph... "); if(bTiming) SW.start();
		P.resize(N, N);
		// get neighboring information
		SymmNeighGraph.resize(N);
		CreateSymmConnectivityGraph(A_OL2, SymmNeighGraph);
		if(bTiming) UG_LOG("took " << SW.ms() << " ms");


	#ifdef UG_PARALLEL
		// 4. coloring in parallel
		//-------------------------------------
		color_process_graph();
		receive_coarsening_from_processes_with_lower_color();
	#endif

		rating.calculate_unassigned();

		// 5. do coarsening
		//-------------------------
		if(0)
		{
			// get possible parent nodes
			UG_LOG("\ncreating possible parent list... "); if(bTiming) SW.start();
			calculate_all_possible_parent_pairs();
			if(bTiming) UG_LOG("took " << SW.ms() << " ms");

			// calculate ratings (not precalculateable because of coarse/uninterpolateable)
			UG_LOG(std::endl << "calculate ratings... "); if(bTiming) SW.start();
			GetRatings(possible_parents, rating, heap);
			if(bTiming) UG_LOG("took " << SW.ms() << " ms");

			heap.print();

			// do coarsening
			UG_LOG(std::endl << "coarsening... "); if(bTiming) SW.start();
			precalculate_coarsening();
			if(bTiming) UG_LOG("took " << SW.ms() << " ms.");
		}
		else
		{
			UG_LOG(std::endl << "other coarsening... "); if(bTiming) SW.start();
			on_demand_coarsening();
			if(bTiming) UG_LOG("took " << SW.ms() << " ms.");
		}

		// [ debug output
		IF_DEBUG(LIB_ALG_AMG, 3)
		{
			UG_LOG("Coarse nodes:\n");
			for(size_t i=0; i<N; i++)
			{
				UG_LOG(i << " rating: " << rating[i]);
				if(rating[i].is_coarse())
					UG_LOG(", new index = " << rating.newIndex[i]);
				UG_LOG("\n");
			}
		}
		write_debug_matrix_markers();
		// ]

		// 6. do aggressive coarsening
		//-----------------------------------------------
		// todo: check if this is ok to be off in parallel

		UG_LOG(std::endl << "second coarsening... "); if(bTiming) SW.start();
		if(m_famg.get_aggressive_coarsening() == true)
			get_aggressive_coarsening_interpolation();
		else
			set_uninterpolateable_as_coarse();
		if(bTiming) UG_LOG("took " << SW.ms() << " ms.");

		#ifdef UG_PARALLEL
		// 7. send coarsening data to processes with higher color
		//----------------------------------------------------------
			send_coarsening_data_to_processes_with_higher_color();
		#endif


		// resize P, output statistics & debug stuff

		UG_LOG(std::endl << N - rating.get_unassigned() << " nodes assigned, " << rating.get_nr_of_coarse() << " coarse, "
				<< N - rating.get_unassigned() - rating.get_nr_of_coarse() << " fine, " << rating.get_unassigned() << " unassigned.");

		// set P to real size
		P.resize(A.num_rows(), rating.get_nr_of_coarse());
		P.defragment();
	#ifdef UG_PARALLEL
		P.set_storage_type(PST_CONSISTENT);
	#endif

		// create parentIndex
		UG_LOG(std::endl << "create parentIndex... "); if(bTiming) SW.start();
		create_parentIndex();
		if(bTiming) UG_LOG("took " << SW.ms() << " ms.");

		//UG_LOG("parentIndex level " << level << "\n")
		//for(size_t i=0; i<rating.get_nr_of_coarse(); i++) { UG_ASSERT(amghelper.parentIndex[level+1][i] != -1, i << " == -1???"); UG_LOG(i << " = " << amghelper.parentIndex[level+1][i] << "\n"); }


		IF_DEBUG(LIB_ALG_AMG, 4)
		{
			UG_LOG("\n\nP:\n");
			P.p();
			UG_LOG("\n\n");
		}



		if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");
		UG_DLOG(LIB_ALG_AMG, 1, std::endl << rating.get_nr_of_coarse() << " coarse, " << N - rating.get_unassigned() - rating.get_nr_of_coarse() << " fine.");

		// 8. construct restriction R = I_{h->2h}
		//-----------------------------------------

		UG_DLOG(LIB_ALG_AMG, 1, std::endl << "restriction... ");	if(bTiming) SW.start();

		UG_LOG(std::endl << "restriction... "); if(bTiming) SW.start();
		// construct restriction R = I_{h -> 2h}
		R.set_as_transpose_of(P);
		// R is already defragmented
		if(bTiming) UG_LOG("took " << SW.ms() << " ms.");

	#ifdef UG_PARALLEL
		R.set_storage_type(PST_CONSISTENT);
	#endif

		if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

	#ifdef FAMG_PRINT_R
		IF_DEBUG(LIB_ALG_AMG, 4)
		{
			UG_LOG(std::endl << "Restriction level " << level << std::endl);
			R.print();
		}
	#endif

		// 9. create Galerkin product AH = R A P
		//-----------------------------------------

		UG_LOG("\ngalerkin product... "); if(bTiming) SW.start();

		// AH = R A P
		CreateAsMultiplyOf(AH, R, A, P);

		if(bTiming) UG_LOG("took " << SW.ms() << " ms");

		// finalize
		if(bTiming) { UG_LOG(std::endl << "Finalizing.."); SW.start(); }
		AH.defragment();
		if(bTiming) UG_LOG(" took " << SW.ms() << " ms\n");

	#ifdef UG_PARALLEL
		AH.set_storage_type(PST_ADDITIVE);
	#endif


	#ifdef FAMG_PRINT_AH
		IF_DEBUG(LIB_ALG_AMG, 3)
		{
			UG_LOG("AH level " << level << std::endl);
			AH.print();
		}
	#endif

		// 9. create Galerkin product AH = R A P
		//-----------------------------------------

		write_debug_matrices();

		// todo: remove dynamic cast, change big_testvector to parallel
		for(size_t i=0; i<m_testvectors.size(); i++)
			CalculateNextTestvector(R, m_testvectors[i]);

		PrintLayout(A_OL2.get_communicator(), nextLevelMasterLayout, nextLevelSlaveLayout);

	}

private:
	void write_debug_matrices();
	void write_debug_matrix_markers();
	template<typename TMatrix>
	void write_debug_matrix(TMatrix &mat, size_t fromlevel, size_t tolevel, const char *name);


	void on_demand_coarsening();
	void precalculate_coarsening();
	void calculate_all_possible_parent_pairs();
};




template<>
void famg<CPUAlgebra>::c_create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
		prolongation_matrix_type &P, size_t level)
{
	FAMG_PROFILE_FUNC();
	stdvector< vector_type > testvectors;
	stdvector<double> omega;

	if(m_testvectors.size() == 0 && m_vVectorWriters.size() == 0)
	{
		testvectors.resize(1);
		omega.resize(1);
		omega[0] = 1.0;
	}
	else
	{
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

			vector_type &vec = testvectors[index];

			std::fstream file((m_writeMatrixPath + "testvector" + ToString(i) + ".values").c_str(),
					std::ios::out);
			for(size_t j=0; j<vec.size(); j++)
				file << j << " " << (vec[j]) << std::endl;

			WriteVectorToConnectionViewer((m_writeMatrixPath + "testvector" + ToString(i) +".mat").c_str(),
					vec, m_amghelper.positions, 2); 

		}
	}

	UG_ASSERT(testvectors.size() > 0, "we need at least one testvector.");

	// testvectors will be altered by FAMGLevelCalculator
	//UG_SET_DEBUG_LEVELS(4);

	FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type> dummy(*this, AH, R, A, P, level, testvectors, omega);
	dummy.do_calculation();

	//UG_SET_DEBUG_LEVELS(0);
}


} // namespace ug

#ifdef UG_PARALLEL
#include "famg_parallel_coarsening_impl.h"
#endif

#include "famg_debug_impl.h"
#include "famg_on_demand_coarsening_impl.h"
#include "famg_precalculate_coarsening_impl.h"

#endif //  __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__
