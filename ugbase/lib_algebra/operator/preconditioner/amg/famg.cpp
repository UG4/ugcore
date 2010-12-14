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

const double omega = 0.6;
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

int famg_log_level = 0;


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

//namespace ug{


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

		calculator.GetAllNeighborsInterpolation(i, P, rating, newIndex, iNrOfCoarse, unassigned);
	}
}


template<typename matrix_type>
void GetPossibleParentPairs(const matrix_type &A, std::vector<std::vector<neighborstruct2> > &possible_neighbors,
		famg_nodes &rating)
{
	possible_neighbors.clear();
	possible_neighbors.resize(A.num_rows());

	UG_LOG(std::scientific);

	FAMGInterpolationCalculator<matrix_type> calculator(A);

	for(size_t i=0; i<A.num_rows(); i++)
		calculator.GetPossibleParentPairs(i, possible_neighbors[i], rating[i]);
}

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

template<typename matrix_type, typename neighborstruct>
void FAMGCoarsen(cgraph &SymmNeighGraph, matrix_type &P, std::vector<std::vector<neighborstruct> > &possible_neighbors,
		famg_nodes &rating, maxheap<famg_nodeinfo> heap, std::vector<int> &newIndex,
		size_t &iNrOfCoarse, size_t &unassigned)
{
	famg_log_level = 4;

	while(heap.height() != 0)
	{
		size_t i = heap.remove_max();
		neighborstruct2 &n = possible_neighbors[i][0];

		FAMG_LOG(2, "\n\n\nSelect next node...\n");
		FAMG_LOG(2, "node " << i << " has rating " << rating[i] << ". now gets fine. parents: ");
		for(size_t j=0; j < n.parents.size(); j++)
			FAMG_LOG(2, n.parents[j].from << " ");

		FAMG_LOG(2, "\nUpdate Neighbors of " << i << "\n");
		// node i gets fine. update neighbors.
		rating[i].set_fine();
		unassigned--;
		UpdateNeighbors(SymmNeighGraph, i, possible_neighbors, rating, heap);

		FAMG_LOG(2, "Set coarse parents:\n");
		// get parent pair, set as coarse (if not already done), update neighbors.

		for(size_t j=0; j < n.parents.size(); j++)
		{
			size_t node = n.parents[j].from;

			if(rating[node].is_coarse()) { FAMG_LOG(2, "\nnode " << node << " is already coarse\n"); }
			else { FAMG_LOG(2, "\nnode " << node << " has ratin " << rating[node] << ". now gets coarse.\nUpdate Neighbors of " << node << "\n"); }

			if(!rating[node].is_coarse())
			{
				heap.remove(node);
				rating[node].set_coarse();
				unassigned--;
				newIndex[node] = iNrOfCoarse++;
				P(node, newIndex[node]) = 1.0;
				UpdateNeighbors(SymmNeighGraph, node, possible_neighbors, rating, heap);
			}
			P(i, newIndex[node]) = n.parents[j].value;
		}
	}

	famg_log_level = 0;

}

template<typename matrix_type>
void FAMGCoarsen2(std::vector<famg_nodeinfo> rating, matrix_type &P, std::vector<int> &newIndex, size_t &iNrOfCoarse)
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
// createFAMGLevel:
//-------------------------
/**
 * create FAMG matrix R, P, and AH = R A P
 * \param AH
 * \param R
 * \param A
 * \param P
 * \param level
 */

void c_create_AMG_level(SparseMatrix<double> &AH, SparseMatrix<double> &R, const SparseMatrix<double> &A,
		SparseMatrix<double> &P, cAMG_helper &amghelper, int level)
{

	bool bTiming=true;
	stopwatch SW, SWwhole; SWwhole.start();

	size_t N = A.num_rows();
	FAMG_LOG(0, "Creating level " << level << ". (" << N << " nodes)" << std::endl << std::fixed);

	std::vector<std::vector<neighborstruct2> > possible_neighbors;
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


	// get possible parent nodes
	FAMG_LOG(2, "creating possible parent list... "); if(bTiming) SW.start();

	GetPossibleParentPairs(A, possible_neighbors, rating);

	if(bTiming) FAMG_LOG(1, "took " << SW.ms() << " ms");


	// calculate ratings (not precalculateable because of coarse/uninterpolateable)
	FAMG_LOG(2, std::endl << "calculate ratings... "); if(bTiming) SW.start();

	GetRatings(possible_neighbors, rating, heap);

	if(bTiming) FAMG_LOG(2, "took " << SW.ms() << " ms");


	// do coarsening
	FAMG_LOG(2, std::endl << "coarsening... "); if(bTiming) SW.start();

	heap.print();
	size_t unassigned = heap.height();

	P.create(N,N);

	// get neighboring information
	cgraph SymmNeighGraph(N);
	CreateSymmConnectivityGraph(A, SymmNeighGraph);

	FAMGCoarsen(SymmNeighGraph, P, possible_neighbors, rating, heap, newIndex, iNrOfCoarse, unassigned);

	IF_FAMG_LOG(2)
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


	FAMG_LOG(2, "\n");

	if(bTiming) FAMG_LOG(2, "took " << SW.ms() << " ms.");
	FAMG_LOG(2, std::endl << N - unassigned << " nodes assigned, " << iNrOfCoarse << " coarse, " << N - unassigned - iNrOfCoarse << " fine, " << unassigned << " unassigned.\n");

	FAMG_LOG(2, std::endl << "second coarsening... "); if(bTiming) SW.start();

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
	UG_LOG("parentIndex level " << level << "\n")
	for(size_t i=0; i<iNrOfCoarse; i++) { UG_ASSERT(amghelper.parentIndex[level+1][i] != -1, i << " == -1???"); UG_LOG(i << " = " << amghelper.parentIndex[level+1][i] << "\n"); }



	FAMG_LOG(2, "\n\nP:\n")
	P.p();
	FAMG_LOG(2, "\n\n");



	if(bTiming) FAMG_LOG(2, "took " << SW.ms() << " ms");
	FAMG_LOG(2, std::endl << iNrOfCoarse << " coarse, " << N - unassigned - iNrOfCoarse << " fine.");

	// construct restriction R = I_{h->2h}
	/////////////////////////////////////////

	FAMG_LOG(2, std::endl << "restriction... ");	if(bTiming) SW.start();

	// construct restriction R = I_{h -> 2h}
	R.create_as_transpose_of(P);
	// R is already finalized

/*#ifdef UG_PARALLEL
	R.set_storage_type(PST_CONSISTENT);
#endif*/

	if(bTiming) FAMG_LOG(2, "took " << SW.ms() << " ms");

#ifdef FAMG_PRINT_R
	IF_FAMG_LOG(3)
	{
		UG_LOG(std::endl << "Restriction level " << level << std::endl);
		R.print();
	}
#endif

	// create Galerkin product
	/////////////////////////////////////////

	FAMG_LOG(2, "\ngalerkin product... "); if(bTiming) SW.start();

	// AH = R A P
	CreateAsMultiplyOf(AH, R, A, P);

	if(bTiming) FAMG_LOG(2, "took " << SW.ms() << " ms");

	// finalize
	if(bTiming) { FAMG_LOG(2, std::endl << "Finalizing.."); SW.start(); }
	AH.finalize();
	if(bTiming) FAMG_LOG(2, " took " << SW.ms() << " ms");

#ifdef FAMG_PRINT_AH
	IF_FAMG_LOG(3)
	{
		UG_LOG("AH level " << level << std::endl);
		AH.print();
	}
#endif

	// finish
	/////////////////////////////////////////

	int nnz = AH.total_num_connections();
	FAMG_LOG(2, std::endl << "AH: nnz: " << nnz << " Density: " <<
			double(nnz)/(double(AH.num_rows())*double(AH.num_rows()))*100.0 << "% nnz/n: " << nnz/(double)AH.num_rows() << std::endl);

	FAMG_LOG(2, "Coarsening rate: " << (100.0*AH.num_rows())/(N) << "%" << std::endl);

	FAMG_LOG(2, " level took " << SWwhole.ms() << " ms" << std::endl);


#ifdef AMG_WRITE_MATRICES_PATH
	if(A.num_rows() < AMG_WRITE_MATRICES_MAX)
	{
		FAMG_LOG(2, "write matrices");
		AMGWriteToFile(A, level, level, (string(AMG_WRITE_MATRICES_PATH) + "A" + ToString(level) + ".mat").c_str(), amghelper);
		FAMG_LOG(2, ".");
		fstream f((string(AMG_WRITE_MATRICES_PATH) + "A" + ToString(level) + ".mat").c_str(), ios::out | ios::app);
		f << "c " << AMG_WRITE_MATRICES_PATH << "fine" << level << ".marks\n";
		f << "c " << AMG_WRITE_MATRICES_PATH << "coarse" << level << ".marks\n";
		f << "c " << AMG_WRITE_MATRICES_PATH << "other" << level << ".marks\n";
		f << "c " << AMG_WRITE_MATRICES_PATH << "dirichlet" << level << ".marks\n";

		AMGWriteToFile(P, level+1, level, (string(AMG_WRITE_MATRICES_PATH) + "Pp" + ToString(level) + ".mat").c_str(), amghelper);
		fstream f2((string(AMG_WRITE_MATRICES_PATH) + "Pp" + ToString(level) + ".mat").c_str(), ios::out | ios::app);
		f2 << "c " << AMG_WRITE_MATRICES_PATH << "fine" << level << ".marks\n";
		f2 << "c " << AMG_WRITE_MATRICES_PATH << "coarse" << level << ".marks\n";
		f2 << "c " << AMG_WRITE_MATRICES_PATH << "other" << level << ".marks\n";
		f2 << "c " << AMG_WRITE_MATRICES_PATH << "dirichlet" << level << ".marks\n";

		FAMG_LOG(2, ".");
		AMGWriteToFile(R, level, level+1, (string(AMG_WRITE_MATRICES_PATH) + "Rr" + ToString(level) + ".mat").c_str(), amghelper);
		fstream f3((string(AMG_WRITE_MATRICES_PATH) + "Rr" + ToString(level) + ".mat").c_str(), ios::out | ios::app);
		f3 << "c " << AMG_WRITE_MATRICES_PATH << "fine" << level << ".marks\n";
		f3 << "c " << AMG_WRITE_MATRICES_PATH << "coarse" << level << ".marks\n";
		f3 << "c " << AMG_WRITE_MATRICES_PATH << "other" << level << ".marks\n";
		f3 << "c " << AMG_WRITE_MATRICES_PATH << "dirichlet" << level << ".marks\n";

		FAMG_LOG(2, ".");
		AMGWriteToFile(AH, level+1, level+1, (string(AMG_WRITE_MATRICES_PATH) + "A" + ToString(level+1) + ".mat").c_str(), amghelper);
		FAMG_LOG(2, ". done.\n");
	}
#endif

	FAMG_LOG(0, "\n");
}


} // namespace ug

#endif //  __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__
