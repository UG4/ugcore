/**
 * \file famg.cpp
 *
 * \author Martin Rupp
 *
 * \date 7.12.2010
 *
 * implementation file for famg
 *
 * Goethe-Center for Scientific Computing 2010.
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

const double omega = 0.0;
const double damping = 0.5;
const double delta = 0.5;


#define AMG_WRITE_MATRICES_PATH "/Users/mrupp/matrices/AMG_"
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

#define IF_FAMG_LOG(level) if(level <= famg_log_level)
#define FAMG_LOG(level, msg) { if(level <= famg_log_level) { UG_LOG(msg); } }

namespace ug{

int famg_log_level = 1;


template<typename value_type>
struct s_interpolation
{
	size_t from;
	value_type value;
};

struct neighborstruct2
{
	void print()
	{
		FAMG_LOG(2, parents.size() << " parents, F = " << F << ": ");
		for(size_t i=0; i<parents.size(); i++)
			FAMG_LOG(2, (i>0 ? "," : "") << "p" << i << ": [" << parents[i].from << "] -> " << parents[i].value);
		FAMG_LOG(2, std::endl);
	}
	FixedArray1<s_interpolation<double>, 2> parents;
	double F;
};

struct neighborstruct_var
{
	void print()
	{
		FAMG_LOG(2, parents.size() << " parents, F = " << F << ": ");
		for(size_t i=0; i<parents.size(); i++)
			FAMG_LOG(2, (i>0 ? "," : "") << "p" << i << ": [" << parents[i].from << "] -> " << parents[i].value);
		FAMG_LOG(2, std::endl);
	}
	std::vector<s_interpolation<double> > parents;
	double F;
};

}
//#define GRAPH_WITH_LOCAL_INVERSE

#include "famg_testvectors.h"
#include "famg_interpolation_calculator.h"
#include "famg_nodeinfo.h"
#include "other_famg_coarsening.h"

//namespace ug{


// GetAggressiveCoarseningInterpolation:
//---------------------------------------
/**
 * \param A				the matrix A
 * \param P				the to-write interpolation matrix P
 * \param rating		rating to determine if a node is uninterpolateable
 * \param newIndex		newIndex[i] = j means node i on the current level has index j on the next coarser
 * \param iNrOfCoarse	number of coarse nodes so far
 * \param unassigned 	number of still unassigned nodes
 */
template<typename matrix_type>
void GetAggressiveCoarseningInterpolation(const matrix_type &A, matrix_type &P,
		famg_nodes &rating, std::vector<int> &newIndex,
				size_t &iNrOfCoarse, size_t &unassigned)
{
	FAMGInterpolationCalculator<matrix_type> calculator(A);

	for(size_t i=0; i<A.num_rows(); i++)
	{
		UG_ASSERT(rating[i].is_valid_rating() == false, "node " << i << " has valid rating (neither coarse nor fine but interpolateable), but is not in heap anymore???");
		if(rating[i].is_uninterpolateable() == false) continue;

		calculator.get_all_neighbors_interpolation(i, P, rating, newIndex, iNrOfCoarse, unassigned);
	}
}


// GetPossibleParentPairs:
//--------------------------
/**
 * \param A						the matrix A
 * \param possible_parents		list of possible interpolating parents for each node
 * \param rating				where the coarsening-ratings for the nodes are stored
 */
template<typename matrix_type>
void GetPossibleParentPairs(const matrix_type &A, std::vector<std::vector<neighborstruct2> > &possible_parents,
		famg_nodes &rating)
{
	possible_parents.clear();
	possible_parents.resize(A.num_rows());

	UG_LOG(std::scientific);

	FAMGInterpolationCalculator<matrix_type> calculator(A);

	for(size_t i=0; i<A.num_rows(); i++)
		calculator.get_possible_parent_pairs(i, possible_parents[i], rating[i]);
}

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
	IF_FAMG_LOG(3)
	{
		UG_LOG("\nSymmNeighGraph:\n");
		SymmNeighGraph.print();
		UG_LOG("\n\n");
	}
#endif
}

// FAMGCoarsen:
//-------------------------
/**
 * performs the FAMG coarsening process with a list of precalculated possible parent nodes
 * \param SymmNeighGraph		used to determine which neighbors' rating needs to be updated
 * \param P						prolongation
 * \param possible_parents		list of possible interpolating parents for each node
 * \param rating				where the coarsening-ratings for the nodes are stored
 * \param heap					heap used for sorting of ratings
 * \param newIndex				newIndex[i] = j means node i on the current level has index j on the next coarser
 * \param iNrOfCoarse			for each coarse node created, this gets increased
 * \param unassigned			for each node assigned coarse or fine, this get decreased
 */
template<typename matrix_type, typename neighborstruct>
void FAMGCoarsen(cgraph &SymmNeighGraph, matrix_type &P, std::vector<std::vector<neighborstruct> > &possible_parents,
		famg_nodes &rating, maxheap<famg_nodeinfo> heap, std::vector<int> &newIndex,
		size_t &iNrOfCoarse, size_t &unassigned)
{

	while(heap.height() != 0)
	{
		// get node i with best rating
		size_t i = heap.remove_max();
		neighborstruct2 &n = possible_parents[i][0];

		FAMG_LOG(2, "\n\n\nSelect next node...\n");
		FAMG_LOG(2, "node " << i << " has rating " << rating[i] << ". now gets fine. parents: ");
		for(size_t j=0; j < n.parents.size(); j++)	FAMG_LOG(2, n.parents[j].from << " ");
		FAMG_LOG(2, "\nUpdate Neighbors of " << i << "\n");

		// node i then gets fine, parent nodes get coarse, and neighbors of those updated.

		// node i gets fine. update neighbors.
		rating[i].set_fine();
		unassigned--;
		UpdateNeighbors(SymmNeighGraph, i, possible_parents, rating, heap);

		FAMG_LOG(2, "Set coarse parents:\n");
		// get parent pair, set as coarse (if not already done), update neighbors.

		for(size_t j=0; j < n.parents.size(); j++)
		{
			size_t node = n.parents[j].from;

			if(rating[node].is_coarse()) { FAMG_LOG(2, "\nnode " << node << " is already coarse\n"); }
			else { FAMG_LOG(2, "\nnode " << node << " has rating " << rating[node] << ". now gets coarse.\nUpdate Neighbors of " << node << "\n"); }

			if(!rating[node].is_coarse())
			{
				heap.remove(node);
				rating[node].set_coarse();
				unassigned--;
				newIndex[node] = iNrOfCoarse++;
				P(node, newIndex[node]) = 1.0;
				UpdateNeighbors(SymmNeighGraph, node, possible_parents, rating, heap);
			}
			P(i, newIndex[node]) = n.parents[j].value;
		}
	}

}

template<typename matrix_type>
void FAMGCoarsen2(famg_nodes &rating, matrix_type &P, std::vector<int> &newIndex, size_t &iNrOfCoarse)
{
	for(size_t i=0; i<rating.size(); i++)
	{
		UG_ASSERT(rating[i].is_valid_rating() == false, "node " << i << " has valid rating (neither coarse nor fine but interpolateable), but is not in heap anymore???");
		if(rating[i].is_uninterpolateable() == false) continue;

		rating[i].set_coarse();
		newIndex[i] = iNrOfCoarse++;
		P(i, newIndex[i]) = 1.0;
	}
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// c_createAMGLevel:
//-------------------------
/**
 * create FAMG matrix R, P, and AH = R A P
 * \param AH
 * \param R
 * \param A
 * \param P
 * \param amghelepr
 * \param level
 */
void c_create_AMG_level(SparseMatrix<double> &AH, SparseMatrix<double> &R, const SparseMatrix<double> &A,
		SparseMatrix<double> &P, cAMG_helper &amghelper, int level)
{

	bool bTiming=true;
	stopwatch SW, SWwhole; SWwhole.start();

	size_t N = A.num_rows();
	FAMG_LOG(0, "Creating level " << level << ". (" << N << " nodes)" << std::endl << std::fixed);

	std::vector<std::vector<neighborstruct2> > possible_parents;
	std::vector<int> newIndex;

	famg_nodes rating(N);

	newIndex.resize(N);
	for(size_t i=0; i<N; i++)
		newIndex[i] = -1;
	size_t iNrOfCoarse=0;

	maxheap<famg_nodeinfo> heap(N, &rating[0]);

#ifdef FAMG_PRINT_A
	IF_FAMG_LOG(3)
	{
		FAMG_LOG("\bThe Matrix A:\n\n");
		A.print();
		FAMG_LOG("\n\n");
	}
#endif
	FAMG_LOG(1, "Create P, SymmNeighGraph... "); if(bTiming) SW.start();
	P.create(N,N);
	// get neighboring information
	cgraph SymmNeighGraph(N);
	CreateSymmConnectivityGraph(A, SymmNeighGraph);
	if(bTiming) FAMG_LOG(1, "took " << SW.ms() << " ms");

	size_t unassigned = N;
	if(0)
	{
		// get possible parent nodes
		FAMG_LOG(1, "creating possible parent list... "); if(bTiming) SW.start();
		GetPossibleParentPairs(A, possible_parents, rating);
		if(bTiming) FAMG_LOG(1, "took " << SW.ms() << " ms");


		// calculate ratings (not precalculateable because of coarse/uninterpolateable)
		FAMG_LOG(2, std::endl << "calculate ratings... "); if(bTiming) SW.start();
		GetRatings(possible_parents, rating, heap);
		if(bTiming) FAMG_LOG(1, "took " << SW.ms() << " ms");

		//rating[40].rating = 1;
		//heap.update(40);
		heap.print();
		unassigned = heap.height();

		// do coarsening
		FAMG_LOG(1, std::endl << "coarsening... "); if(bTiming) SW.start();
		FAMGCoarsen(SymmNeighGraph, P, possible_parents, rating, heap, newIndex, iNrOfCoarse, unassigned);
		if(bTiming) FAMG_LOG(1, "took " << SW.ms() << " ms.");
	}
	else
	{
		FAMG_LOG(1, std::endl << "other coarsening... "); if(bTiming) SW.start();
		other_coarsening(A, SymmNeighGraph, P, possible_parents, rating, heap, newIndex, iNrOfCoarse, unassigned);
		if(bTiming) FAMG_LOG(1, "took " << SW.ms() << " ms.");
	}


	IF_FAMG_LOG(4)
	{
		UG_LOG("Coarse nodes:\n");
		for(size_t i=0; i<N; i++)
		{
			//if(newIndex[i] != -1)
			UG_LOG(i << " rating: " << rating[i]);
			if(newIndex[i] != -1)
				UG_LOG(", new index = " << newIndex[i]);
			UG_LOG("\n");
		}
	}


	{
	fstream ffine((string(AMG_WRITE_MATRICES_PATH) + "fine" + ToString(level) + ".marks").c_str(), std::ios::out);
	fstream fcoarse((string(AMG_WRITE_MATRICES_PATH) + "coarse" + ToString(level) + ".marks").c_str(), std::ios::out);
	fstream fother((string(AMG_WRITE_MATRICES_PATH) + "other" + ToString(level) + ".marks").c_str(), std::ios::out);
	fstream fdirichlet((string(AMG_WRITE_MATRICES_PATH) + "dirichlet" + ToString(level) + ".marks").c_str(), std::ios::out);
	for(size_t i=0; i<N; i++)
	{
		int o = amghelper.GetOriginalIndex(level, i);
		if(rating[i].is_fine()) ffine << o << "\n";
		else if(rating[i].is_coarse()) fcoarse << o << "\n";
		else if(rating[i].is_dirichlet()) fdirichlet << o << "\n";
		else fother << o << "\n";
	}
	}


	FAMG_LOG(1, std::endl << N - unassigned << " nodes assigned, " << iNrOfCoarse << " coarse, "
			<< N - unassigned - iNrOfCoarse << " fine, " << unassigned << " unassigned.");

	FAMG_LOG(1, std::endl << "second coarsening... "); if(bTiming) SW.start();

	//FAMGCoarsen2(rating, P, newIndex, iNrOfCoarse);
	GetAggressiveCoarseningInterpolation(A, P, rating, newIndex, iNrOfCoarse, unassigned);


	P.resize(N, iNrOfCoarse);

	// create parentIndex
	amghelper.parentIndex[level+1] = new int[iNrOfCoarse];
	for(size_t i=0; i<iNrOfCoarse; i++) amghelper.parentIndex[level+1][i] = -1;
	for(size_t i=0; i<N; i++)
	{
		if(newIndex[i] != -1)
			amghelper.parentIndex[level+1][ newIndex[i] ] = i;

	}
	//UG_LOG("parentIndex level " << level << "\n")
	//for(size_t i=0; i<iNrOfCoarse; i++) { UG_ASSERT(amghelper.parentIndex[level+1][i] != -1, i << " == -1???"); UG_LOG(i << " = " << amghelper.parentIndex[level+1][i] << "\n"); }


	IF_FAMG_LOG(4)
	{
		FAMG_LOG(4, "\n\nP:\n")
		P.p();
		FAMG_LOG(4, "\n\n");
	}



	if(bTiming) FAMG_LOG(1, "took " << SW.ms() << " ms");
	FAMG_LOG(1, std::endl << iNrOfCoarse << " coarse, " << N - unassigned - iNrOfCoarse << " fine.");

	// construct restriction R = I_{h->2h}
	/////////////////////////////////////////

	FAMG_LOG(1, std::endl << "restriction... ");	if(bTiming) SW.start();

	// construct restriction R = I_{h -> 2h}
	R.create_as_transpose_of(P);
	// R is already finalized

/*#ifdef UG_PARALLEL
	R.set_storage_type(PST_CONSISTENT);
#endif*/

	if(bTiming) FAMG_LOG(1, "took " << SW.ms() << " ms");

#ifdef FAMG_PRINT_R
	IF_FAMG_LOG(3)
	{
		UG_LOG(std::endl << "Restriction level " << level << std::endl);
		R.print();
	}
#endif

	// create Galerkin product
	/////////////////////////////////////////

	FAMG_LOG(1, "\ngalerkin product... "); if(bTiming) SW.start();

	// AH = R A P
	CreateAsMultiplyOf(AH, R, A, P);

	if(bTiming) FAMG_LOG(1, "took " << SW.ms() << " ms");

	// finalize
	if(bTiming) { FAMG_LOG(1, std::endl << "Finalizing.."); SW.start(); }
	//AH.finalize();
	if(bTiming) FAMG_LOG(1, " took " << SW.ms() << " ms");

#ifdef FAMG_PRINT_AH
	IF_FAMG_LOG(3)
	{
		UG_LOG("AH level " << level << std::endl);
		AH.print();
	}
#endif

	// finish
	/////////////////////////////////////////

	FAMG_LOG(1, " level took " << SWwhole.ms() << " ms" << std::endl);


#ifdef AMG_WRITE_MATRICES_PATH
	if(A.num_rows() < AMG_WRITE_MATRICES_MAX)
	{
		FAMG_LOG(1, "write matrices");
		AMGWriteToFile(A, level, level, (string(AMG_WRITE_MATRICES_PATH) + "A" + ToString(level) + ".mat").c_str(), amghelper);
		fstream f((string(AMG_WRITE_MATRICES_PATH) + "A" + ToString(level) + ".mat").c_str(), ios::out | ios::app);
		f << "c " << AMG_WRITE_MATRICES_PATH << "fine" << level << ".marks\n";
		f << "c " << AMG_WRITE_MATRICES_PATH << "coarse" << level << ".marks\n";
		f << "c " << AMG_WRITE_MATRICES_PATH << "other" << level << ".marks\n";
		f << "c " << AMG_WRITE_MATRICES_PATH << "dirichlet" << level << ".marks\n";
		FAMG_LOG(1, ".");

		AMGWriteToFile(P, level+1, level, (string(AMG_WRITE_MATRICES_PATH) + "Pp" + ToString(level) + ".mat").c_str(), amghelper);
		fstream f2((string(AMG_WRITE_MATRICES_PATH) + "Pp" + ToString(level) + ".mat").c_str(), ios::out | ios::app);
		f2 << "c " << AMG_WRITE_MATRICES_PATH << "fine" << level << ".marks\n";
		f2 << "c " << AMG_WRITE_MATRICES_PATH << "coarse" << level << ".marks\n";
		f2 << "c " << AMG_WRITE_MATRICES_PATH << "other" << level << ".marks\n";
		f2 << "c " << AMG_WRITE_MATRICES_PATH << "dirichlet" << level << ".marks\n";
		FAMG_LOG(1, ".");

		AMGWriteToFile(R, level, level+1, (string(AMG_WRITE_MATRICES_PATH) + "Rr" + ToString(level) + ".mat").c_str(), amghelper);
		fstream f3((string(AMG_WRITE_MATRICES_PATH) + "Rr" + ToString(level) + ".mat").c_str(), ios::out | ios::app);
		f3 << "c " << AMG_WRITE_MATRICES_PATH << "fine" << level << ".marks\n";
		f3 << "c " << AMG_WRITE_MATRICES_PATH << "coarse" << level << ".marks\n";
		f3 << "c " << AMG_WRITE_MATRICES_PATH << "other" << level << ".marks\n";
		f3 << "c " << AMG_WRITE_MATRICES_PATH << "dirichlet" << level << ".marks\n";
		FAMG_LOG(1, ".");

		AMGWriteToFile(AH, level+1, level+1, (string(AMG_WRITE_MATRICES_PATH) + "A" + ToString(level+1) + ".mat").c_str(), amghelper);
		FAMG_LOG(1, ". done.\n");
	}
#endif
}


} // namespace ug

#endif //  __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__
