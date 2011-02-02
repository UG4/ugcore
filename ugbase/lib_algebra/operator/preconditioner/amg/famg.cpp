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
#include "famg.h"


#include "lib_algebra/parallelization/parallel_matrix_overlap_impl.h"

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

namespace ug{

/*#ifdef UG_PARALLEL
typedef ParallelMatrix<SparseMatrix<double> > matrix_type;
#else
typedef SparseMatrix<double> matrix_type;
#endif*/


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
		UG_LOG(parents.size() << " parents, F = " << F << ": ");
		for(size_t i=0; i<parents.size(); i++)
			UG_LOG((i>0 ? "," : "") << "p" << i << ": [" << GetOriginalIndex(parents[i].from) << "] -> " << parents[i].value);
		UG_LOG(std::endl);
	}
	FixedArray1<s_interpolation<double>, 2> parents;
	double F;
};

struct neighborstruct_var
{
	void print()
	{
		UG_LOG(parents.size() << " parents, F = " << F << ": ");
		for(size_t i=0; i<parents.size(); i++)
			UG_LOG((i>0 ? "," : "") << "p" << i << ": [" << GetOriginalIndex(parents[i].from) << "] -> " << parents[i].value);
		UG_LOG(std::endl);
	}
	stdvector<s_interpolation<double> > parents;
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
		famg_nodes &rating, stdvector<int> &newIndex,
				size_t &iNrOfCoarse, size_t &unassigned,
				FAMGInterpolationCalculator<matrix_type> &calculator)
{
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
void GetPossibleParentPairs(const matrix_type &A, stdvector<stdvector<neighborstruct2> > &possible_parents,
		famg_nodes &rating, FAMGInterpolationCalculator<matrix_type> &calculator)
{
	possible_parents.clear();
	possible_parents.resize(A.num_rows());

	UG_LOG(std::scientific);

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

	IF_DEBUG(LIB_ALG_AMG, 3)
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
void FAMGCoarsen(cgraph &SymmNeighGraph, matrix_type &P, stdvector<stdvector<neighborstruct> > &possible_parents,
		famg_nodes &rating, maxheap<famg_nodeinfo> heap, stdvector<int> &newIndex,
		size_t &iNrOfCoarse, size_t &unassigned)
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
		unassigned--;
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
void SetUninterpolateableAsCoarse(famg_nodes &rating, matrix_type &P, stdvector<int> &newIndex, size_t &iNrOfCoarse)
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

///	algebra blocks are static.
template <> struct block_traits<MathVector<3> >
{
	enum{
		is_static = 1
	};
};


int ColorProcessorGraph(pcl::ParallelCommunicator<IndexLayout> &com, std::set<int> &pids)
{
	// 1. send all neighbors our number of neighbors
	size_t myDegree = pids.size();
	typedef std::set<int>::iterator setiterator;
	std::map<int, size_t> othersDegree;

	UG_DLOG(LIB_ALG_AMG, 2, "sending degree data to ");
	for(setiterator iter = pids.begin(); iter != pids.end(); ++iter)
	{
		int pid = *iter;
		com.send_raw(pid, &myDegree, sizeof(size_t), true);
		com.receive_raw(pid, &othersDegree[pid], sizeof(size_t));
		UG_DLOG(LIB_ALG_AMG, 2, pid << " ");
	}
	UG_DLOG(LIB_ALG_AMG, 2, "\n");

	com.communicate();

	int myPID = pcl::GetProcRank();
	UG_DLOG(LIB_ALG_AMG, 2, "process " <<  myPID << " here. I have degree " << myDegree << ".\n");
	for(setiterator iter = pids.begin(); iter != pids.end(); ++iter)
	{
		int pid = *iter;
		UG_DLOG(LIB_ALG_AMG, 2, " neighbor pid " << pid << " has degree " << othersDegree[pid] << ".\n");
	}

	stdvector<int> colors(myDegree, -1);

	size_t processesWithHigherWeight=0, processesWithLowerWeight=0;

	// 2. for all neighbors which have a higher weight than this process, issue a color receive
	UG_DLOG(LIB_ALG_AMG, 2, "issue receive from... ");
	for(setiterator iter = pids.begin(); iter != pids.end(); ++iter)
	{
		int pid2 = *iter; size_t degree2 = othersDegree[pid2];
		if(myDegree < degree2 || (myDegree == degree2 && myPID < pid2))
		{
			com.receive_raw(pid2, &colors[processesWithHigherWeight], sizeof(int));
			processesWithHigherWeight++;
			UG_DLOG(LIB_ALG_AMG, 2, pid2 << " ");
		}
		else
			processesWithLowerWeight++;
	}
	UG_DLOG(LIB_ALG_AMG, 2, "\n");

	int myColor = 0;

	// 3. communicate & get a free color (only if actually got color from neighbors)
	if(processesWithHigherWeight > 0)
	{
		UG_DLOG(LIB_ALG_AMG, 2, "communicate...");
		com.communicate();
		UG_DLOG(LIB_ALG_AMG, 2, "done.\n");

		UG_DLOG(LIB_ALG_AMG, 2, "received colors from " << processesWithHigherWeight << " processes!\n");
		size_t i=0;
		for(setiterator iter = pids.begin(); iter != pids.end(); ++iter)
		{
			int pid2 = *iter; size_t degree2 = othersDegree[pid2];
			if(myDegree < degree2 || (myDegree == degree2 && myPID < pid2))
			{
				if(colors[i] == -1)
				UG_DLOG(LIB_ALG_AMG, 2, "pid " << pid2 << " did not send color?");
				UG_DLOG(LIB_ALG_AMG, 2, " neighbor pid " << pid2 << " has color " << colors[i] << ".\n");
				i++;
			}
		}
		colors.resize(processesWithHigherWeight);

		sort(colors.begin(), colors.end());
		for(; myColor < (int)colors.size(); myColor++)
		{
			if(myColor != colors[myColor])
				break;
		}
	}

	UG_DLOG(LIB_ALG_AMG, 2, "GOT COLOR! My color is " << myColor << "\n");

	// 4. for all neighbors which have a lower weight than this process, issue a color send

	if(processesWithLowerWeight > 0)
	{
		UG_DLOG(LIB_ALG_AMG, 2, "Sending color to " << processesWithLowerWeight << " processes: ")
		for(setiterator iter = pids.begin(); iter != pids.end(); ++iter)
		{
			int pid2 = *iter; size_t degree2 = othersDegree[pid2];
			if(myDegree > degree2 || (myDegree == degree2 && myPID > pid2 ))
			{
				// bug when fixedsize=false???
				com.send_raw(pid2, &myColor, sizeof(int), true);
				UG_DLOG(LIB_ALG_AMG, 2, pid2 << " ");
			}
		}
		UG_DLOG(LIB_ALG_AMG, 2, "communicate...");
		com.communicate();
		UG_DLOG(LIB_ALG_AMG, 2, "done.\n");
	}

	UG_DLOG(LIB_ALG_AMG, 2, "coloring done.\n");

	return myColor;
}

template<>
void famg<CPUAlgebra>::c_create_AMG_level(matrix_type &AH, SparseMatrix<double> &R, const matrix_type &A,
		SparseMatrix<double> &P, size_t level)
{
	//UG_SET_DEBUG_LEVELS(4);
	UG_LOG("Creating level " << level << ". (" << A.num_rows() << " nodes)" << std::endl << std::fixed);
//	bool bTiming=true;
	stopwatch SW, SWwhole; SWwhole.start();



	if(level == 0)
	{
		matrix_type A_OL2;
		std::vector<IndexLayout> masterLayouts, slaveLayouts;
		IndexLayout totalMasterLayout, totalSlaveLayout;

		GenerateOverlap(A, A_OL2, totalMasterLayout, totalSlaveLayout, masterLayouts, slaveLayouts, 2);

		std::vector<MathVector<3> > vec2(&amghelper.positions[0], &amghelper.positions[amghelper.size]);
		vec2.resize(A_OL2.num_rows());

		/*UG_DLOG(LIB_ALG_AMG, 0, "positions:\n");
		for(size_t i=0; i<vec2.size(); ++i)
			UG_DLOG(LIB_ALG_MATRIX, 0, i << " = " << vec2[i] << "\n");*/


	//	copy all ids from master to slave interfaces
		ComPol_VecCopy<std::vector<MathVector<3> > >	copyPol(&vec2);

		pcl::ParallelCommunicator<IndexLayout> &communicator = A_OL2.get_communicator();
		communicator.send_data(totalMasterLayout, copyPol);
		communicator.receive_data(totalSlaveLayout, copyPol);
		communicator.communicate();

		UG_DLOG(LIB_ALG_AMG, 0, "positions after copy:\n");
		for(size_t i=0; i<vec2.size(); ++i)
			UG_DLOG(LIB_ALG_MATRIX, 0, i << " = " << vec2[i] << "\n");

		std::stringstream ss; ss << "A_OL2_" << pcl::GetProcRank() << ".mat";
		WriteMatrixToConnectionViewer(ss.str().c_str(), A_OL2, &vec2[0], 2);


		for(size_t i=0; i<slaveLayouts.size(); i++)
		{
			std::set<int> s;
			for(IndexLayout::iterator iter = slaveLayouts[i].begin(); iter != slaveLayouts[i].end(); ++iter)
			{
				IndexLayout::Interface &interface = slaveLayouts[i].interface(iter);
				for(IndexLayout::Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
					s.insert(interface.get_element(iter2));
			}

			UG_LOG("OL Level " << i << ": ");
			for(std::set<int>::iterator iter = s.begin(); iter != s.end(); ++iter)
			{
				UG_LOG(*iter << " ");
			}
			UG_LOG("\n");
		}

		std::set<int> pids;
		for(IndexLayout::iterator iter = totalMasterLayout.begin(); iter != totalMasterLayout.end(); ++iter)
			pids.insert(totalMasterLayout.proc_id(iter));
		for(IndexLayout::iterator iter = totalSlaveLayout.begin(); iter != totalSlaveLayout.end(); ++iter)
			pids.insert(totalSlaveLayout.proc_id(iter));

		int myColor = ColorProcessorGraph(communicator, pids);
	}

}
#if 0

	FAMGInterpolationCalculator<matrix_type> calculator(A, m_delta, m_theta, m_dDampingForSmootherInInterpolationCalculation);

	UG_LOG("calculating testvector... ");
	if(bTiming) SW.start();
	if(big_testvector.size() == 0)
	{
		big_testvector.resize(A.num_rows());
		for(size_t i=0; i<A.num_rows(); i++)
			if(A.is_isolated(i) && m_bTestvectorZeroAtDirichlet)
				big_testvector[i] = 0.0;
			else
				big_testvector[i] = 1.0;

		//big_testvector = 1.0;
	}

	{
		Vector<double> d; d.resize(A.num_rows());
		for(size_t jj=0; jj < m_iTestvectorDamps; jj++)
		{
			MatMult(d, 1.0, A, big_testvector);
			for(size_t i=0; i<A.num_rows(); i++)
				big_testvector[i] = big_testvector[i] - 0.6*d[i]/A(i,i);
		}
	}

	if(bTiming) UG_LOG("took " << SW.ms() << " ms");

	pAmghelper = &amghelper;
	currentlevel = level;

	size_t N = A.num_rows();

	stdvector<stdvector<neighborstruct2> > possible_parents;
	stdvector<int> newIndex;

	famg_nodes rating(N);

#ifdef UG_PARALLEL
	// mark slaves


#endif


	newIndex.resize(N);	for(size_t i=0; i<N; i++)	newIndex[i] = -1;
	size_t iNrOfCoarse=0;

	maxheap<famg_nodeinfo> heap(N, &rating[0]);

	IF_DEBUG(LIB_ALG_AMG, 3)
	{
		UG_LOG("\bThe Matrix A:\n\n");
		A.print();
		UG_LOG("\n\n");
	}

	UG_LOG("\nCreate P, SymmNeighGraph... "); if(bTiming) SW.start();
	P.create(N,N);
	// get neighboring information
	cgraph SymmNeighGraph(N);
	CreateSymmConnectivityGraph(A, SymmNeighGraph);
	if(bTiming) UG_LOG("took " << SW.ms() << " ms");

	size_t unassigned = N;
	if(0)
	{
		// get possible parent nodes
		UG_LOG("\ncreating possible parent list... "); if(bTiming) SW.start();
		GetPossibleParentPairs(A, possible_parents, rating, calculator);
		if(bTiming) UG_LOG("took " << SW.ms() << " ms");


		// calculate ratings (not precalculateable because of coarse/uninterpolateable)
		UG_LOG(std::endl << "calculate ratings... "); if(bTiming) SW.start();
		GetRatings(possible_parents, rating, heap);
		if(bTiming) UG_LOG("took " << SW.ms() << " ms");

		//rating[40].rating = 1;
		//heap.update(40);
		heap.print();
		unassigned = heap.height();

		// do coarsening
		UG_LOG(std::endl << "coarsening... "); if(bTiming) SW.start();
		FAMGCoarsen(SymmNeighGraph, P, possible_parents, rating, heap, newIndex, iNrOfCoarse, unassigned);
		if(bTiming) UG_LOG("took " << SW.ms() << " ms.");
	}
	else
	{
		UG_LOG(std::endl << "other coarsening... "); if(bTiming) SW.start();
		other_coarsening(A, SymmNeighGraph, P, possible_parents, rating, heap, newIndex, iNrOfCoarse, unassigned, calculator);
		if(bTiming) UG_LOG("took " << SW.ms() << " ms.");
	}


	IF_DEBUG(LIB_ALG_AMG, 4)
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


	if(m_writeMatrices)
	{
		fstream ffine((m_writeMatrixPath + "AMG_fine" + ToString(level) + ".marks").c_str(), std::ios::out);
		fstream fcoarse((m_writeMatrixPath + "AMG_coarse" + ToString(level) + ".marks").c_str(), std::ios::out);
		fstream fother((m_writeMatrixPath + "AMG_other" + ToString(level) + ".marks").c_str(), std::ios::out);
		fstream fdirichlet((m_writeMatrixPath + "AMG_dirichlet" + ToString(level) + ".marks").c_str(), std::ios::out);
		for(size_t i=0; i<N; i++)
		{
			int o = amghelper.GetOriginalIndex(level, i);
			if(rating[i].is_fine()) ffine << o << "\n";
			else if(rating[i].is_coarse()) fcoarse << o << "\n";
			else if(rating[i].is_dirichlet()) fdirichlet << o << "\n";
			else fother << o << "\n";
		}
	}


	UG_LOG(std::endl << N - unassigned << " nodes assigned, " << iNrOfCoarse << " coarse, "
			<< N - unassigned - iNrOfCoarse << " fine, " << unassigned << " unassigned.");

	UG_LOG(std::endl << "second coarsening... "); if(bTiming) SW.start();


	if(m_bAggressiveCoarsening)
		GetAggressiveCoarseningInterpolation(A, P, rating, newIndex, iNrOfCoarse, unassigned, calculator);
	else
		SetUninterpolateableAsCoarse(rating, P, newIndex, iNrOfCoarse);

	if(bTiming) UG_LOG("took " << SW.ms() << " ms.");

	P.resize(N, iNrOfCoarse);
	P.finalize();

	// create parentIndex
	UG_LOG(std::endl << "create parentIndex... "); if(bTiming) SW.start();
	stdvector<stdvector<int> > &parentIndex = *amghelper.parentIndex;
	parentIndex.resize(level+2);
	parentIndex[level+1].resize(iNrOfCoarse);
	for(size_t i=0; i<iNrOfCoarse; i++) parentIndex[level+1][i] = -1;
	for(size_t i=0; i<N; i++) {
		if(newIndex[i] != -1) parentIndex[level+1][ newIndex[i] ] = i;
	}
	if(bTiming) UG_LOG("took " << SW.ms() << " ms.");

	//UG_LOG("parentIndex level " << level << "\n")
	//for(size_t i=0; i<iNrOfCoarse; i++) { UG_ASSERT(amghelper.parentIndex[level+1][i] != -1, i << " == -1???"); UG_LOG(i << " = " << amghelper.parentIndex[level+1][i] << "\n"); }


	IF_DEBUG(LIB_ALG_AMG, 4)
	{
		UG_LOG("\n\nP:\n");
		P.p();
		UG_LOG("\n\n");
	}



	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");
	UG_DLOG(LIB_ALG_AMG, 1, std::endl << iNrOfCoarse << " coarse, " << N - unassigned - iNrOfCoarse << " fine.");

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

	if(m_writeMatrices && A.num_rows() < AMG_WRITE_MATRICES_MAX)
	{
		UG_DLOG(LIB_ALG_AMG, 1, "write matrices");
		AMGWriteToFile(A, level, level, (m_writeMatrixPath + "AMG_A" + ToString(level) + ".mat").c_str(), amghelper);
		fstream f((m_writeMatrixPath + "AMG_A" + ToString(level) + ".mat").c_str(), ios::out | ios::app);
		f << "c " << m_writeMatrixPath << "AMG_fine" << level << ".marks\n";
		f << "c " << m_writeMatrixPath << "AMG_coarse" << level << ".marks\n";
		f << "c " << m_writeMatrixPath << "AMG_other" << level << ".marks\n";
		f << "c " << m_writeMatrixPath << "AMG_dirichlet" << level << ".marks\n";
		f << "v " << m_writeMatrixPath << "AMG_d" << level << ".values\n";
		UG_DLOG(LIB_ALG_AMG, 1, ".");

		AMGWriteToFile(P, level+1, level, (m_writeMatrixPath + "AMG_Pp" + ToString(level) + ".mat").c_str(), amghelper);
		fstream f2((m_writeMatrixPath + "AMG_Pp" + ToString(level) + ".mat").c_str(), ios::out | ios::app);
		f2 << "c " << m_writeMatrixPath << "AMG_fine" << level << ".marks\n";
		f2 << "c " << m_writeMatrixPath << "AMG_coarse" << level << ".marks\n";
		f2 << "c " << m_writeMatrixPath << "AMG_other" << level << ".marks\n";
		f2 << "c " << m_writeMatrixPath << "AMG_dirichlet" << level << ".marks\n";
		UG_DLOG(LIB_ALG_AMG, 1, ".");

		AMGWriteToFile(R, level, level+1, (m_writeMatrixPath + "AMG_Rr" + ToString(level) + ".mat").c_str(), amghelper);
		fstream f3((m_writeMatrixPath + "AMG_Rr" + ToString(level) + ".mat").c_str(), ios::out | ios::app);
		f3 << "c " << m_writeMatrixPath << "AMG_fine" << level << ".marks\n";
		f3 << "c " << m_writeMatrixPath << "AMG_coarse" << level << ".marks\n";
		f3 << "c " << m_writeMatrixPath << "AMG_other" << level << ".marks\n";
		f3 << "c " << m_writeMatrixPath << "AMG_dirichlet" << level << ".marks\n";
		UG_DLOG(LIB_ALG_AMG, 1, ".");

		AMGWriteToFile(AH, level+1, level+1, (m_writeMatrixPath + "AMG_A" + ToString(level+1) + ".mat").c_str(), amghelper);
		UG_DLOG(LIB_ALG_AMG, 1, ". done.\n");
	}


	Vector<double> t;
	t.resize(R.num_rows());

	MatMult(t, 1.0, R, big_testvector);
	big_testvector.resize(R.num_rows());
	big_testvector = t;
}
#endif

} // namespace ug

#endif //  __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__


