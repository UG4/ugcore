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

#include "amg_nodeinfo.h"
#include "stopwatch.h"


using namespace std;
namespace ug{

//#define GRAPH_WITH_LOCAL_INVERSE




#define AMG_WRITE_MATRICES_PATH "/Users/mrupp/matrices/AMG_"
#define AMG_WRITE_MATRICES_MAX (200*200)

#define AMG_WRITE_GRAPH

#if 0

#define LATE_COARSE_SOLVER // do coarsening down to 10 nodes.
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

template<typename T> inline double amg_diag_value(const T &d) { return d.norm(); }
template<typename T> inline double amg_offdiag_value(const T &d) { return -d.norm(); }


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
template<typename Matrix_type>
void CreateStrongConnectionGraph(const Matrix_type &A, cgraph &graph, double theta)
{
	graph.resize(A.num_rows());

	for(size_t i=0; i< A.num_rows(); i++)
	{
		if(A[i].is_isolated())
			continue;

		double dmax = 0;

		for(typename Matrix_type::cRowIterator conn = A.beginRow(i); !conn.isEnd(); ++conn)
		{
			if((*conn).iIndex == i) continue; // skip diag
			if((*conn).dValue != 0.0 && amg_offdiag_value((*conn).dValue) < dmax)
				dmax = amg_offdiag_value((*conn).dValue);
		}

		for(typename Matrix_type::cRowIterator conn = A.beginRow(i); !conn.isEnd(); ++conn)
		{
			if((*conn).iIndex == i) continue; // skip diag
			if( amg_offdiag_value((*conn).dValue) < theta * dmax)
				graph.set_connection(i, (*conn).iIndex);
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
template<typename Matrix_type, typename Vector_type>
void amg<Matrix_type, Vector_type>::createAMGLevel(Matrix_type &AH, SparseMatrix<double> &R, const Matrix_type &A, SparseMatrix<double> &P, int level)
{
	amg_nodeinfo *nodes = new amg_nodeinfo[A.num_rows()];

	bool bTiming=true;
	cout << "Creating level " << level << ". (" << A.num_rows() << " nodes)" << endl;
	stopwatch SW;
	stopwatch SWwhole; SWwhole.start();

	// amg_nodeinfo: infos zu den einzelnen knoten von A, verwaltung der ratings etc

	// todo: check for isolated condition

	maxheap<amg_nodeinfo> PQ;

	//std::vector<bool coarse(A.num_rows());
	int unassigned = A.num_rows();

	int *newIndex = new int[A.num_rows()];
	for(size_t i=0; i<A.num_rows(); ++i) newIndex[i] = -1;

	int *posInConnections = new int[A.num_rows()];
	for(size_t i=0; i<A.num_rows(); ++i) posInConnections[i] = -1;


	cout << "A.totalNrOfConnections = " << A.total_num_connections() << endl;

	cgraph graphS(A.num_rows());
	cgraph graphST;

	// build graph
	/////////////////////////////////////////

	cout << "building graph... "; cout.flush();
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


	CreateMeasureOfImportancePQ(graphS, graphST, PQ, unassigned, nodes);

	if(bTiming) SW.printTimeDiff();

#ifdef AMG_PRINT_COARSEN_RATINGS
	 for(size_t i=0; i<A.num_rows(); i++)
		 cout << i << " [" << GetOriginalIndex(level, i) << "] " << nodes[i] << endl;
#endif
	//PQ.print();

	// Coarsen
	/////////////////////////////////////////
	cout << "coarsening... "; cout.flush();

	if(bTiming) SW.start();
	int iNrOfCoarse = 0;
	Coarsen(graphST, PQ, newIndex, unassigned, iNrOfCoarse, nodes);
	PreventFFConnections(graphS, graphST, nodes, newIndex, iNrOfCoarse);

	if(bTiming) { SW.printTimeDiff();}

	// agressive Coarsening
	/////////////////////////////////////////

	if(aggressiveCoarsening && level == 0)
	{
		// build graph 2
		//------------------

		cgraph graph2(A.num_rows());

		cout << "building graph2... "; cout.flush();
		if(bTiming) SW.start();

		for(size_t i=0; i < A.num_rows(); i++) newIndex[i] = -1;

		cgraph graphAC;
		CreateAggressiveCoarseningGraph(graphST, graphAC, nodes, aggressiveCoarseningNrOfPaths, posInConnections);
		CreateMeasureOfImportanceAggressiveCoarseningPQ(graphAC, PQ, unassigned, iNrOfCoarse, newIndex, nodes);
		if(bTiming) SW.printTimeDiff();

		/*for(int i=0; i<A.num_rows(); i++)
		 {
		 cout << i << " (" << GetPosForIndexAtLevel(i, level).x << " " << GetPosForIndexAtLevel(i, level).y << ") " ;
		 nodes[i].print();
		 }//*/
		//PQ.print();

		// coarsen 2
		//------------------

		if(unassigned == 0)
			cout << "skipping coarsening2: no unassigned nodes." << endl;
		else
		{
			cout << "coarsening2... "; cout.flush();

			if(bTiming) SW.start();
			Coarsen(graphAC, PQ, newIndex, unassigned, iNrOfCoarse, nodes);
			//PreventFFConnections(graphS, graphST, nodes, newIndex, iNrOfCoarse);
			if(bTiming) { SW.printTimeDiff();}

		}
	}

	// create vectors for AMG multigrid
	/////////////////////////////////////////

	vec1[level+1] = new Vector_type (iNrOfCoarse);
	vec2[level+1] = new Vector_type (iNrOfCoarse);
	cout << "created vec1 on level" << level +1 << endl;

	// set size for variable sized blockvectors
	for(size_t i=0; i<A.num_rows(); i++)
		if(nodes[i].isCoarse())
		{
			int rows = GetRows((*A.beginRow(i)).dValue);
			UG_ASSERT(newIndex[i] >= 0, "");
			SetSize((*vec1[level+1])[newIndex[i]], rows);
			SetSize((*vec2[level+1])[newIndex[i]], rows);
		}

	// set parentindex for debugging
//#ifdef FLEXAMG
	parentIndex[level+1] = new int[iNrOfCoarse];
	for(size_t i=0; i<A.num_rows(); i++)
		if(nodes[i].isCoarse())
			parentIndex[level+1][ newIndex[i] ] = i;
//#endif

	/*for(size_t i=0; i<A.num_rows(); i++)
	 {
	 cout << i << " (" << GetPosForIndexAtLevel(i, level).x << " " << GetPosForIndexAtLevel(i, level).y << ") ";
	 nodes[i].print();
	 }//*/

#ifdef AMG_PRINT_COARSENING
	printCoarsening(level);
#endif

	// construct prolongation P = I_{2h->h}
	/////////////////////////////////////////

	cout << "prolongation... "; cout.flush();

	unassigned = 0;

#if FLEXAMG
	P.fromlevel = level+1;
	P.tolevel = level;
#endif

	if(bTiming) SW.start();
	CreateRugeStuebenProlongation(P, A, newIndex, iNrOfCoarse, unassigned, nodes, theta);
	//UG_ASSERT(unassigned == 0 || (aggressiveCoarsening && level == 0),
		//	"no aggressive coarsening, but indirect interpolation of " << unassigned " nodes ?");
	if(unassigned > 0)
		CreateIndirectProlongation(P, A, newIndex, unassigned, nodes, posInConnections, theta);

	if(bTiming) SW.printTimeDiff();

#ifdef AMG_PRINT_P
	cout << "Prolongation level " << level << endl;
	P.print();
#endif

	// construct prolongation R = I_{h->2h}
	/////////////////////////////////////////

	cout << "restriction... "; cout.flush();
	if(bTiming) SW.start();
	// construct restriction R = I_{h -> 2h}
	R.create_as_transpose_of(P); // already finished
	//R.name = "AMG:R";
	//R.print("R");
	if(bTiming) SW.printTimeDiff();

#ifdef AMG_PRINT_R
	cout << "Restriction level " << level << endl;
	R.print();
#endif

	// create Galerkin product
	/////////////////////////////////////////

	cout << "galerkin product... "; cout.flush();
	if(bTiming) SW.start();

	// AH = R A P
	CreateAsMultiplyOf(AH, R, A, P);

	if(bTiming) SW.printTimeDiff();

	if(bTiming) SW.start();
	AH.finalize();
	if(bTiming) { cout << "finalizing..."; SW.printTimeDiff(); }

#ifdef AMG_PRINT_AH
	cout << "AH level " << level << endl;
	AH.print();
#endif

	// finish
	/////////////////////////////////////////

	int nnz = AH.total_num_connections();
	cout << "AH: nnz: " << nnz << " Density: " << double(nnz)/(double(AH.num_rows())*double(AH.num_rows()))*100.0 << "% nnz/n: " << nnz/(double)AH.num_rows() << endl;
	cout << "Coarsening rate: " << (100.0*AH.num_rows())/(A.num_rows()) << "%" << endl;

	cout << " level "; SWwhole.printTimeDiff();  cout << endl;
	cout << endl;
	cout.flush();


#ifdef AMG_WRITE_MATRICES_PATH
	if(this->A[0]->num_rows() < AMG_WRITE_MATRICES_MAX)
	{
		cout << "write matrices";
		AMGWriteToFile(P, level+1, level, (string(AMG_WRITE_MATRICES_PATH) + "P" + ToString(level) + ".mat").c_str(), amghelper); cout << "."; cout.flush();
		AMGWriteToFile(R, level, level+1, (string(AMG_WRITE_MATRICES_PATH) + "R" + ToString(level) + ".mat").c_str(), amghelper); cout << "."; cout.flush();
		AMGWriteToFile(AH, level+1, level+1, (string(AMG_WRITE_MATRICES_PATH) + "A" + ToString(level+1) + ".mat").c_str(), amghelper); cout << "."; cout.flush();
		cout << " done." << endl;
	}
#endif

	//P.print();

	cout << AH << endl;
	cout << endl;//*/

	delete[] posInConnections;
	delete[] newIndex;
	delete[] nodes;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// amg<Matrix_type, Vector_type>::init
//----------------
//! creates MG Hierachy for with Matrix_type A and temporary vectors for higher levels
//! @param A	matrix A.
template<typename Matrix_type, typename Vector_type>
bool amg<Matrix_type, Vector_type>::init(const Matrix_type& A_)
{
	amghelper.positions = &dbg_positions[0];
	amghelper.size = A_.num_rows();
	amghelper.parentIndex = parentIndex;
	amghelper.dimension = dbg_dimension;

	cout << "Starting AMG Setup." << endl << endl;
	stopwatch SW;
	const Matrix_type *pA = &A_;
	A[0] = const_cast<Matrix_type*> (pA);

#ifdef FLEXAMG
	A[0]->fromlevel = 0;
	A[0]->tolevel = 0;
#endif


#ifdef AMG_WRITE_MATRICES_PATH
	if(A[0]->num_rows() < AMG_WRITE_MATRICES_MAX)
	{
		cout << "write matrix A...";
		AMGWriteToFile(*A[0], 0, 0, (string(AMG_WRITE_MATRICES_PATH) + "A0.mat").c_str(), amghelper);
		cout << "done." << endl; cout.flush();
	}
#endif

	int i=0;
	while(i< max_levels-1)
	{

		double L = A[i]->num_rows();
#ifndef LATE_COARSE_SOLVER
		//if(L < 100 || A[i]->total_num_connections()/(L*L) > 0.5)	break; // abbruch falls density > 50%
		if(L < 1000)	break; // abbruch falls density > 50%
#else
		if(L < 10)	break;
#endif
		//smoother[i].init(*A[i]);

		A[i+1] = new Matrix_type();
		createAMGLevel(*A[i+1], R[i], *A[i], P[i], i);

		vec3[i] = new Vector_type (A[i]->num_rows());
		i++;
	}

	int nrOfUnknowns = block_vector_traits< typename Vector_type::entry_type >::nrOfUnknowns;
	cout << "Creating level " << i << " (" << A[i]->num_rows() << " nodes, total " << A[i]->num_rows()*nrOfUnknowns << " unknowns)" << endl << "Using Direct Solver on Matrix "
	<< A[i]->num_rows()*nrOfUnknowns << "x" << A[i]->num_rows()*nrOfUnknowns << ". ";

#ifdef FLEXAMG
	stopwatch SW2; SW2.start();
	coarseSolver.create(*A[i]);
	SW2.printTimeDiff();
	cout << endl << endl;
#else
	stopwatch SW2; SW2.start();
	coarseSolver.init(*A[i]);
	SW2.printTimeDiff();
	cout << endl << endl;
#endif

	used_levels = i+1;
	cout << "AMG Setup finished. Used Levels: " << used_levels << ". ";
	SW.printTimeDiff();

	// calc complexities
	double nnzs=0;
	double totallength=0;
	for(int i=0; i<used_levels; i++)
	{
		nnzs += A[i]->total_num_connections();
		totallength += A[i]->num_rows();
	}

	cout << "Operator Complexity: " << nnzs/A[0]->total_num_connections() << " nodes complexity: " << totallength/A[0]->num_rows() << endl << endl;

	return true;
}

//!
//! amg constructor
template<typename Matrix_type, typename Vector_type>
amg<Matrix_type, Vector_type>::amg()
{
	used_levels = 0;
	max_levels = 10;
	aggressiveCoarsening = 0;
	aggressiveCoarseningNrOfPaths = 2; // A2

	nu1 = 2;
	nu2 = 2;
	gamma = 1;

	eps_truncation_of_interpolation = 0.3; // no truncation (or 0.2).

	sigma = 0.3;
	theta = 0.3;

	//FORCE_CREATION { printCoarsening(0,0); }
}

//!
//! amg destructor
template<typename Matrix_type, typename Vector_type>
amg<Matrix_type, Vector_type>::~amg()
{
	for(int i=1; i<used_levels-1; i++)
	{
		delete A[i];
		delete vec1[i];
		delete vec2[i];
		delete vec3[i-1];
		if(parentIndex[i]) delete [] parentIndex[i];
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// get_correction_and_update_defect:
//------------------------------------

template<typename Matrix_type, typename Vector_type>
bool amg<Matrix_type, Vector_type>::get_correction_and_update_defect(Vector_type &d, Vector_type &c, int level)
{
	UG_ASSERT(c.size() == d.size() && c.size() == A[level]->num_rows(),
			"c.size = " << c.size() << ", d.size = " << d.size() << ", A.size = " << A[level]->num_rows() << ": not matching");

	const Matrix_type &Ah = *(A[level]);

	if(level == used_levels-1)
	{
		coarseSolver.apply(d, c);
		d -= Ah*c;
		//Ah.matmul_minus(d, c);
		return 0.1e-14;
	}

	Vector_type &corr = *vec3[level];

	// presmooth
	for(int i=0; i < nu1; i++)
	{
		// calc c = B^{-1} b
		diag_step(Ah, corr, d, 0.66);
		c += corr;
		d -= Ah*corr;
	}

	Vector_type &cH = *vec1[level+1];
	Vector_type &dH = *vec2[level+1];

	// restrict defect
	dH = R[level]*d;
	//R[level].apply(dH, d);

	cH = 0.0;

	// apply lmgc on coarser nodes
	if(level+1 == used_levels-1)
		get_correction_and_update_defect(dH, cH, level+1);
	else
		for(int i=0; i<gamma; i++)
			get_correction_and_update_defect(dH, cH, level+1);

	//interpolate correction
	corr = P[level]*cH;
	//P[level].apply(corr, cH);

	// add coarse grid correction to level correction
	c += corr;

	//update defect
	d -= Ah*corr;

	// postsmooth
	for(int i=0; i < nu2; i++)
	{
		diag_step(Ah, corr, d, 0.66);
		c += corr;
		d-= Ah*corr;
		//A.matmul_minus(d, c);
	}

	return true;
}

} // namespace ug

#endif //  __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_IMPL_H__
