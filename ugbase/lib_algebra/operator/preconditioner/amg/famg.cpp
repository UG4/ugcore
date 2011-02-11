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
 * for test purposes, functions are here in a cpp file.
 */


#ifndef __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__
#define __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__

//#include "sparsematrix_util.h"

#include "ug.h"
#include "amg_debug_helper.h"

#include "stopwatch.h"
#include "common/assert.h"
#include "maxheap.h"
#include "famg.h"
#include <set>

#ifdef UG_PARALLEL
std::string GetFilename(std::string path, std::string name, std::string extension)
{
	return path + name + "_" + ToString(pcl::GetProcRank()) + extension;
}
#else
string GetFilename(string path, string name, string extension)
{
	return path + name + "_0" + extension;
}
#endif

ug::Vector<double> big_testvector;
ug::cAMG_helper *pAmghelper;
int currentlevel=-1;


size_t GetOriginalIndex(size_t i)
{
	UG_ASSERT(currentlevel != -1, "");
	return pAmghelper->GetOriginalIndex(currentlevel, i);
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


int ColorProcessorGraph(pcl::ParallelCommunicator<IndexLayout> &com, std::set<int> &pids,
		std::vector<int> &processesWithLowerColor,
		std::vector<int> &processesWithHigherColor);
/*#ifdef UG_PARALLEL
typedef ParallelMatrix<SparseMatrix<double> > matrix_type;
#else
typedef SparseMatrix<double> matrix_type;
#endif*/


// CreateSymmConnectivityGraph:
//------------------------------
/**
 * \param A						the matrix A
 * \param SymmNeighGraph		later used to determine which neighbors' rating needs to be updated
 */
template<typename matrix_type>
void CreateSymmConnectivityGraph(const matrix_type &A, cgraph &SymmNeighGraph)
{
	for(size_t i=0; i<A.num_rows(); i++)
		for(typename matrix_type::cRowIterator conn = A.beginRow(i); !conn.isEnd(); ++conn)
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




template<typename matrix_type, typename prolongation_matrix_type>
class FAMGLevelCalculator
{
private:
	famg<CPUAlgebra> m_famg;
	matrix_type &AH;
	const matrix_type &A;
	bool bTiming;


	prolongation_matrix_type &R, &P;
	size_t level;

	FAMGInterpolationCalculator<matrix_type> calculator;
	famg_nodes rating;
	stdvector<stdvector<neighborstruct2> > possible_parents; // list of possible interpolating parents for each node
	maxheap<famg_nodeinfo> heap;  // heap used for sorting of ratings
	cgraph SymmNeighGraph; ///< used to determine which neighbors' rating needs to be updated

#ifdef UG_PARALLEL
	IndexLayout OLReceivingLayout, OLSendingLayout;
	std::vector<int> processesWithLowerColor, processesWithHigherColor;
	int m_myColor;
#endif

#ifdef UG_PARALLEL
	matrix_type A_OL2;
	IndexLayout &nextLevelMasterLayout;
		IndexLayout &nextLevelSlaveLayout;

#else
	matrix_type &A_OL2;
#endif



public:
	FAMGLevelCalculator(famg<CPUAlgebra> &f,
			matrix_type &_AH, prolongation_matrix_type &_R,  const matrix_type &_A,
			prolongation_matrix_type &_P, size_t _level)
	: m_famg(f), AH(_AH), A(_A), R(_R), P(_P), level(_level),
			calculator(A_OL2, m_famg.m_delta, m_famg.m_theta, m_famg.m_dDampingForSmootherInInterpolationCalculation), rating(P)
#ifndef UG_PARALLEL
			, A_OL2(A)
#else
			, nextLevelMasterLayout(*(new IndexLayout)), nextLevelSlaveLayout(*(new IndexLayout)) // TODO: fix that, this is stupid
#endif
	{
	}

private:
	void get_aggressive_coarsening_interpolation()
	{
		for(size_t i=0; i<A.num_rows(); i++)
		{
			UG_ASSERT(rating[i].is_valid_rating() == false, "node " << i << " has valid rating (neither coarse nor fine but interpolateable), but is not in heap anymore???");
			if(rating[i].is_uninterpolateable() == false) continue;

			calculator.get_all_neighbors_interpolation(i, P, rating);
		}
	}


	void calculate_all_possible_parent_pairs(stdvector<stdvector<neighborstruct2> > &possible_parents,
			famg_nodes &rating, FAMGInterpolationCalculator<matrix_type> &calculator)
	{
		possible_parents.clear();
		possible_parents.resize(A.num_rows());

		UG_LOG(std::scientific);

		for(size_t i=0; i<A.num_rows(); i++)
			calculator.get_possible_parent_pairs(i, possible_parents[i], rating);
	}



	void old_coarsen()
	{
		while(heap.height() != 0)
		{
			// get node i with best rating
			size_t i = heap.remove_max();
			neighborstruct2 &n = possible_parents[i][0];

			UG_DLOG(LIB_ALG_AMG, 2, "\n\n\nSelect next node...\n");
			UG_DLOG(LIB_ALG_AMG, 2, "node " << i << " has rating " << rating[i] << ". now gets fine. parents: ");
			for(size_t j=0; j < n.parents.size(); j++)	UG_DLOG(LIB_ALG_AMG, 2, n.parents[j].from << " ");
			UG_DLOG(LIB_ALG_AMG, 2, "\nUpdate Neighbors of " << i << "\n");

			// node i then gets fine, parent nodes get coarse, and neighbors of those updated.

			// node i gets fine. update neighbors.
			rating[i].set_fine();
			UpdateNeighbors(SymmNeighGraph, i, possible_parents, rating, heap);

			UG_DLOG(LIB_ALG_AMG, 2, "Set coarse parents:\n");
			// get parent pair, set as coarse (if not already done), update neighbors.

			for(size_t j=0; j < n.parents.size(); j++)
			{
				size_t node = n.parents[j].from;

				if(rating[node].is_coarse()) { UG_DLOG(LIB_ALG_AMG, 2, "\nnode " << node << " is already coarse\n"); }
				else { UG_DLOG(LIB_ALG_AMG, 2, "\nnode " << node << " has rating " << rating[node] << ". now gets coarse.\nUpdate Neighbors of " << node << "\n"); }

				if(!rating[node].is_coarse())
				{
					heap.remove(node);
					rating.set_coarse(node);
					UpdateNeighbors(SymmNeighGraph, node, possible_parents, rating, heap);
				}
				P(i, rating.newIndex[node]) = n.parents[j].value;
			}
		}

	}

	void set_uninterpolateable_as_coarse()
	{
		for(size_t i=0; i<rating.size(); i++)
		{
			UG_ASSERT(rating[i].is_valid_rating() == false, "node " << i << " has valid rating (neither coarse nor fine but interpolateable), but is not in heap anymore???");
			if(rating[i].is_uninterpolateable() == false) continue;
			rating.set_coarse(i);
		}
	}


#ifdef UG_PARALLEL
	void color_process_graph();
	void receive_coarsening_from_processes_with_lower_color();
	void send_coarsening_data_to_processes_with_higher_color();
	void create_OL2_matrix();
#endif

	void do_stuff()
	{
		UG_LOG("Creating level " << level << ". (" << A.num_rows() << " nodes)" << std::fixed);
		bTiming=true;
		stopwatch SW, SWwhole; SWwhole.start();

		int me = pcl::GetProcRank();
		UG_ASSERT(pcl::GetNumProcesses() < 2, "only for 2, change layouts so that all slaves of a master are informed about a coarse node");



#ifdef UG_PARALLEL
		create_OL2_matrix();
#endif
		size_t N = A_OL2.num_rows();
		rating.create(N);

		UG_LOG("\ncalculating testvector... ");
		if(bTiming) SW.start();

		CalculateTestvector(A_OL2, big_testvector, m_famg.m_bTestvectorZeroAtDirichlet,
				m_famg.m_iTestvectorDamps);

		if(bTiming) UG_LOG("took " << SW.ms() << " ms");


				size_t iNrOfCoarse=0;

		heap.create(N, &rating[0]);

		IF_DEBUG(LIB_ALG_AMG, 3)
		{
			UG_LOG("\bThe Matrix A:\n\n");
			A.print();
			UG_LOG("\n\n");
		}

		UG_LOG("\nCreate P, SymmNeighGraph... "); if(bTiming) SW.start();
		P.create(A.num_rows(),A.num_rows());
		// get neighboring information
		cgraph SymmNeighGraph(N);
		CreateSymmConnectivityGraph(A_OL2, SymmNeighGraph);
		if(bTiming) UG_LOG("took " << SW.ms() << " ms");



	#ifdef UG_PARALLEL
		color_process_graph();
		receive_coarsening_from_processes_with_lower_color();
	#endif

		rating.calculate_unassigned();

		UG_SET_DEBUG_LEVELS(4);
		if(0)
		{
			// get possible parent nodes
			UG_LOG("\ncreating possible parent list... "); if(bTiming) SW.start();
			calculate_all_possible_parent_pairs();
			if(bTiming) UG_LOG("took " << SW.ms() << " ms");


			// calculate ratings (not precalculateable because of coarse/uninterpolateable)
			UG_LOG(std::endl << "calculate ratings... "); if(bTiming) SW.start();
			calculate_all_possible_parent_pairs();
			if(bTiming) UG_LOG("took " << SW.ms() << " ms");

			heap.print();

			// do coarsening
			UG_LOG(std::endl << "coarsening... "); if(bTiming) SW.start();
			old_coarsen();
			if(bTiming) UG_LOG("took " << SW.ms() << " ms.");
		}
		else
		{
			UG_LOG(std::endl << "other coarsening... "); if(bTiming) SW.start();
			other_coarsening();
			if(bTiming) UG_LOG("took " << SW.ms() << " ms.");
		}

	#ifdef UG_PARALLEL
		send_coarsening_data_to_processes_with_higher_color();
	#endif

		IF_DEBUG(LIB_ALG_AMG, 4)
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


		UG_LOG(std::endl << N - rating.unassigned << " nodes assigned, " << rating.iNrOfCoarse << " coarse, "
				<< N - rating.unassigned - rating.iNrOfCoarse << " fine, " << rating.unassigned << " unassigned.");

		UG_LOG(std::endl << "second coarsening... "); if(bTiming) SW.start();


		if(m_famg.m_bAggressiveCoarsening)
			get_aggressive_coarsening_interpolation();
		else
			set_uninterpolateable_as_coarse();

		if(bTiming) UG_LOG("took " << SW.ms() << " ms.");

		P.resize(A.num_rows(), rating.iNrOfCoarse);
		P.finalize();

		// create parentIndex
		UG_LOG(std::endl << "create parentIndex... "); if(bTiming) SW.start();
		stdvector<stdvector<int> > &parentIndex = *m_famg.amghelper.parentIndex;
		parentIndex.resize(level+2);
		parentIndex[level+1].resize(rating.iNrOfCoarse);
		for(size_t i=0; i < rating.iNrOfCoarse; i++) parentIndex[level+1][i] = -1;
		for(size_t i=0; i<N; i++) {
			if(rating.newIndex[i] != -1) parentIndex[level+1][ rating.newIndex[i] ] = i;
		}
		if(bTiming) UG_LOG("took " << SW.ms() << " ms.");

		//UG_LOG("parentIndex level " << level << "\n")
		//for(size_t i=0; i<rating.iNrOfCoarse; i++) { UG_ASSERT(amghelper.parentIndex[level+1][i] != -1, i << " == -1???"); UG_LOG(i << " = " << amghelper.parentIndex[level+1][i] << "\n"); }


		IF_DEBUG(LIB_ALG_AMG, 4)
		{
			UG_LOG("\n\nP:\n");
			P.p();
			UG_LOG("\n\n");
		}



		if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");
		UG_DLOG(LIB_ALG_AMG, 1, std::endl << iNrOfCoarse << " coarse, " << N - rating.unassigned - iNrOfCoarse << " fine.");

		// construct restriction R = I_{h->2h}
		/////////////////////////////////////////

		UG_DLOG(LIB_ALG_AMG, 1, std::endl << "restriction... ");	if(bTiming) SW.start();

		UG_LOG(std::endl << "restriction... "); if(bTiming) SW.start();
		// construct restriction R = I_{h -> 2h}
		R.create_as_transpose_of(P);
		// R is already finalized
		if(bTiming) UG_LOG("took " << SW.ms() << " ms.");

	/*#ifdef UG_PARALLEL
		R.set_storage_type(PST_CONSISTENT);
	#endif*/

		if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");

	#ifdef FAMG_PRINT_R
		IF_DEBUG(LIB_ALG_AMG, 4)
		{
			UG_LOG(std::endl << "Restriction level " << level << std::endl);
			R.print();
		}
	#endif

		// create Galerkin product
		/////////////////////////////////////////

		UG_LOG("\ngalerkin product... "); if(bTiming) SW.start();

		// AH = R A P
		CreateAsMultiplyOf(AH, R, A, P);

		if(bTiming) UG_LOG("took " << SW.ms() << " ms");

		// finalize
		if(bTiming) { UG_LOG(std::endl << "Finalizing.."); SW.start(); }
		AH.finalize();
		if(bTiming) UG_LOG(" took " << SW.ms() << " ms\n");

	#ifdef FAMG_PRINT_AH
		IF_DEBUG(LIB_ALG_AMG, 3)
		{
			UG_LOG("AH level " << level << std::endl);
			AH.print();
		}
	#endif

		// finish
		/////////////////////////////////////////

		write_debug_matrices();


		CalculateNextTestvector(R, big_testvector);

	}

	void write_debug_matrices()
	{
		if(m_famg.m_writeMatrices && A.num_rows() < AMG_WRITE_MATRICES_MAX)
		{
			UG_DLOG(LIB_ALG_AMG, 1, "write matrices");

			write_debug_matrix(A, "AMG_A", level);			UG_DLOG(LIB_ALG_AMG, 1, ".");
			write_debug_matrix(P, "AMG_P", level);			UG_DLOG(LIB_ALG_AMG, 1, ".");
			write_debug_matrix(R, "AMG_R", level);			UG_DLOG(LIB_ALG_AMG, 1, ".");

			AMGWriteToFile(AH, level+1, level+1, GetFilename(m_famg.m_writeMatrixPath, ToString("AMG_A") + ToString(level+1),".mat").c_str(), m_famg.amghelper);

			UG_DLOG(LIB_ALG_AMG, 1, ". done.\n");
		}
	}

	void write_debug_matrix_markers()
	{
		int me = pcl::GetProcRank();
		if(m_famg.m_writeMatrices)
		{
			std::fstream ffine((m_famg.m_writeMatrixPath + "AMG_fine" + ToString(level)
					+ "_" + ToString(me) + ".marks").c_str(), std::ios::out);
			std::fstream fcoarse((m_famg.m_writeMatrixPath + "AMG_coarse" + ToString(level)
					+ "_" + ToString(me) + ".marks").c_str(), std::ios::out);
			std::fstream fother((m_famg.m_writeMatrixPath + "AMG_other" + ToString(level)
					+ "_" + ToString(me) + ".marks").c_str(), std::ios::out);
			std::fstream fdirichlet((m_famg.m_writeMatrixPath + "AMG_dirichlet" + ToString(level)
					+ "_" + ToString(me) + ".marks").c_str(), std::ios::out);
			for(size_t i=0; i < rating.size(); i++)
			{
				int o = m_famg.amghelper.GetOriginalIndex(level, i);
				if(rating[i].is_fine()) ffine << o << "\n";
				else if(rating[i].is_coarse()) fcoarse << o << "\n";
				else if(rating[i].is_dirichlet()) fdirichlet << o << "\n";
				else fother << o << "\n";
			}
		}

	}

	template<typename TMatrix>
	void write_debug_matrix(TMatrix &mat, size_t level, const char *name)
	{
		std::string filename = GetFilename(m_famg.m_writeMatrixPath, ToString(name) + ToString(level),".mat");
		AMGWriteToFile(A, level, level, filename.c_str(), m_famg.amghelper);
		std::fstream f2(filename, std::ios::out | std::ios::app);
		f2 << "c " << GetFilename(m_famg.m_writeMatrixPath, "AMG_fine" + ToString(level), ".marks") << "\n";
		f2 << "c " << GetFilename(m_famg.m_writeMatrixPath, "AMG_coarse" + ToString(level), ".marks") << "\n";
		f2 << "c " << GetFilename(m_famg.m_writeMatrixPath, "AMG_other" + ToString(level), ".marks") << "\n";
		f2 << "c " << GetFilename(m_famg.m_writeMatrixPath, "AMG_dirichlet" + ToString(level), ".marks") << "\n";
	}

	void other_coarsening();
};


template<>
void famg<CPUAlgebra>::c_create_AMG_level(matrix_type &AH, SparseMatrix<double> &R, const matrix_type &A,
		SparseMatrix<double> &P, size_t level)
{
	pAmghelper = &amghelper;
	currentlevel = level;

	FAMGLevelCalculator<matrix_type, SparseMatrix<double> > dummy(*this, AH, R, A, P, level);
}


} // namespace ug

#include "famg_parallel_coarsening_impl.h"
#include "other_famg_coarsening.h"

#endif //  __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__


