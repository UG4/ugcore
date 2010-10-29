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

#define LATE_COARSE_SOLVER // do coarsening down to 10 nodes.


#if 0

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
template<typename matrix_type>
void CreateStrongConnectionGraph(const matrix_type &A, cgraph &graph, double theta)
{
	graph.resize(A.num_rows());

	for(size_t i=0; i< A.num_rows(); i++)
	{
		if(A[i].is_isolated())
			continue;

		double dmax = 0;

		for(typename matrix_type::cRowIterator conn = A.beginRow(i); !conn.isEnd(); ++conn)
		{
			if((*conn).iIndex == i) continue; // skip diag
			if((*conn).dValue != 0.0 && amg_offdiag_value((*conn).dValue) < dmax)
				dmax = amg_offdiag_value((*conn).dValue);
		}

		for(typename matrix_type::cRowIterator conn = A.beginRow(i); !conn.isEnd(); ++conn)
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
template<typename TAlgebra>
void amg<TAlgebra>::create_AMG_level(matrix_type &AH, SparseMatrix<double> &R, const matrix_type &A,
		SparseMatrix<double> &P, int level)
{
	amg_nodeinfo *nodes = new amg_nodeinfo[A.num_rows()];

	bool bTiming=true;
	UG_LOG("Creating level " << level << ". (" << A.num_rows() << " nodes)" << endl << std::fixed);
	stopwatch SW;
	stopwatch SWwhole; SWwhole.start();

	// amg_nodeinfo: infos zu den einzelnen knoten von A, verwaltung der ratings etc

	// todo: check for isolated condition

	nodeinfo_pq_type PQ;

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

	if(bTiming) UG_LOG("took " << SW.ms() << " ms");

#ifdef AMG_PRINT_COARSEN_RATINGS
	 for(size_t i=0; i<A.num_rows(); i++)
		 cout << i << " [" << GetOriginalIndex(level, i) << "] " << nodes[i] << endl;
#endif

	// Coarsen
	/////////////////////////////////////////
	cout << endl << "coarsening... "; cout.flush();

	SW.start();
	int iNrOfCoarse = 0;

	Coarsen(graphST, PQ, newIndex, unassigned, iNrOfCoarse, nodes);

	PreventFFConnections(graphS, graphST, nodes, newIndex, iNrOfCoarse);

	if(bTiming) UG_LOG("took " << SW.ms() << " ms");

	// agressive Coarsening
	/////////////////////////////////////////

	if(aggressiveCoarsening && level == 0)
	{
		// build graph 2
		//------------------

		cgraph graph2(A.num_rows());

		cout << endl << "building graph2... "; cout.flush();
		if(bTiming) SW.start();

		for(size_t i=0; i < A.num_rows(); i++) newIndex[i] = -1;

		cgraph graphAC;
		CreateAggressiveCoarseningGraph(graphST, graphAC, nodes, aggressiveCoarseningNrOfPaths, posInConnections);
		CreateMeasureOfImportanceAggressiveCoarseningPQ(graphAC, PQ, unassigned, iNrOfCoarse, newIndex, nodes);
		if(bTiming) UG_LOG("took " << SW.ms() << " ms");
		// coarsen 2
		//------------------

		if(unassigned == 0)
		{
			UG_LOG(endl << "skipping coarsening2: no unassigned nodes.");
		}
		else
		{
			UG_LOG(endl << "coarsening2... ");

			if(bTiming) SW.start();
			Coarsen(graphAC, PQ, newIndex, unassigned, iNrOfCoarse, nodes);
			//PreventFFConnections(graphS, graphST, nodes, newIndex, iNrOfCoarse);
			if(bTiming) UG_LOG("took " << SW.ms() << " ms");

		}
	}

	// create vectors for AMG multigrid
	/////////////////////////////////////////

	vec3[level] = new vector_type;
	vec3[level]->create(A.num_rows());

	vec1[level+1] = new vector_type;
	vec1[level+1]->create(iNrOfCoarse);
	vec2[level+1] = new vector_type;
	vec2[level+1]->create(iNrOfCoarse);

#ifdef UG_PARALLEL
	// todo: change this for later "right" parallel implementation
	UG_ASSERT(pcl::GetNumProcesses() == 1, "AMG currently only for 1 process");

	//typedef pcl::SingleLevelLayout<pcl::OrderedInterface<size_t, std::vector> > IndexLayout

	vec3[level]->set_communicator(*com[level]);
	// todo: implement "real" ParallelCommunicator
	com[level+1] = new pcl::ParallelCommunicator<IndexLayout>;
	vec1[level+1]->set_communicator(*com[level+1]);
	vec2[level+1]->set_communicator(*com[level+1]);

	vec3[level]->set_storage_type(PST_ADDITIVE);
	vec1[level+1]->set_storage_type(PST_ADDITIVE);
	vec2[level+1]->set_storage_type(PST_ADDITIVE);


	vec3[level]->set_master_layout(pseudoLayout);
	vec3[level]->set_slave_layout(pseudoLayout);
	vec1[level+1]->set_master_layout(pseudoLayout);
	vec1[level+1]->set_slave_layout(pseudoLayout);
	vec2[level+1]->set_master_layout(pseudoLayout);
	vec2[level+1]->set_slave_layout(pseudoLayout);

#endif


	UG_LOG(endl << "created vec3 on level " << level << ", vec1 and vec2 on level" << level +1);

	// todo: set size for variable sized blockvectors
	/*for(size_t i=0; i<A.num_rows(); i++)
		if(nodes[i].isCoarse())
		{
			int rows = GetRows((*A.beginRow(i)).dValue);
			UG_ASSERT(newIndex[i] >= 0, "");
			SetSize((*vec1[level+1])[newIndex[i]], rows);
			SetSize((*vec2[level+1])[newIndex[i]], rows);
		}*/

	// set parentindex for debugging
	parentIndex[level+1] = new int[iNrOfCoarse];
	for(size_t i=0; i<A.num_rows(); i++)
		if(nodes[i].isCoarse())
			parentIndex[level+1][ newIndex[i] ] = i;


#ifdef AMG_PRINT_COARSENING
	printCoarsening(level);
#endif

	// construct prolongation P = I_{2h->h}
	/////////////////////////////////////////

	UG_LOG(endl << "prolongation... ");
	unassigned = 0;

	if(bTiming) SW.start();
	CreateRugeStuebenProlongation(P, A, newIndex, iNrOfCoarse, unassigned, nodes, theta);
	//UG_ASSERT(unassigned == 0 || (aggressiveCoarsening && level == 0),
		//	"no aggressive coarsening, but indirect interpolation of " << unassigned " nodes ?");
	if(unassigned > 0)
		CreateIndirectProlongation(P, A, newIndex, unassigned, nodes, posInConnections, theta);

/*#ifdef UG_PARALLEL
	P.set_storage_type(PST_CONSISTENT);
#endif*/

	if(bTiming) UG_LOG("took " << SW.ms() << " ms");

#ifdef AMG_PRINT_P
	cout << endl << "Prolongation level " << level << endl;
	P.print();
#endif


	// construct restriction R = I_{h->2h}
	/////////////////////////////////////////

	UG_LOG(endl << "restriction... ");
	if(bTiming) SW.start();

	// construct restriction R = I_{h -> 2h}
	R.create_as_transpose_of(P);
	// R is already finalized

/*#ifdef UG_PARALLEL
	R.set_storage_type(PST_CONSISTENT);
#endif*/

	if(bTiming) UG_LOG("took " << SW.ms() << " ms");

#ifdef AMG_PRINT_R
	cout << endl << "Restriction level " << level << endl;
	R.print();
#endif

	// create Galerkin product
	/////////////////////////////////////////

	UG_LOG("\ngalerkin product... ");
	if(bTiming) SW.start();

	// AH = R A P
	CreateAsMultiplyOf(AH, R, A, P);

	if(bTiming) UG_LOG("took " << SW.ms() << " ms");

	// finalize
	if(bTiming) SW.start();
	AH.finalize();
	if(bTiming) UG_LOG(endl << "Finalizing... took " << SW.ms() << " ms");

#ifdef AMG_PRINT_AH
	cout << "AH level " << level << endl;
	AH.print();
#endif

	// finish
	/////////////////////////////////////////

	int nnz = AH.total_num_connections();
	UG_LOG(endl << "AH: nnz: " << nnz << " Density: " <<
			double(nnz)/(double(AH.num_rows())*double(AH.num_rows()))*100.0 << "% nnz/n: " << nnz/(double)AH.num_rows() << endl);

	UG_LOG("Coarsening rate: " << (100.0*AH.num_rows())/(A.num_rows()) << "%" << endl);

	UG_LOG(" level took " << SWwhole.ms() << " ms" << endl);


#ifdef AMG_WRITE_MATRICES_PATH
	if(this->A[0]->num_rows() < AMG_WRITE_MATRICES_MAX)
	{
		UG_LOG("write matrices");
		AMGWriteToFile(P, level+1, level, (string(AMG_WRITE_MATRICES_PATH) + "P" + ToString(level) + ".mat").c_str(), amghelper);
		UG_LOG(".");
		AMGWriteToFile(R, level, level+1, (string(AMG_WRITE_MATRICES_PATH) + "R" + ToString(level) + ".mat").c_str(), amghelper);
		UG_LOG(".");
		AMGWriteToFile(AH, level+1, level+1, (string(AMG_WRITE_MATRICES_PATH) + "A" + ToString(level+1) + ".mat").c_str(), amghelper);
		UG_LOG(". done.\n");
	}
#endif

	UG_LOG("\n");



	delete[] posInConnections;
	delete[] newIndex;
	delete[] nodes;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// amg<matrix_type, vector_type>::init
//----------------
//! creates MG Hierachy for with matrix_type A and temporary vectors for higher levels
//! @param A	matrix A.
template<typename TAlgebra>
bool amg<TAlgebra>::preprocess(matrix_type& mat)
{
	cleanup();
	A[0] = &mat;
	m_bInited = false;
	return true;
}

template<typename TAlgebra>
bool amg<TAlgebra>::init_amg()
{
	if(m_basesolver==NULL)
	{
		UG_LOG("AMG::init_amg(): No base solver selected. Call set_base_solver(basesolver) to set a base solver.\n");
		return false;
	}
	if(m_presmoother==NULL)
	{
		UG_LOG("AMG::init_amg(): No PreSmoother selected. Call set_presmoother(presmoother) to set a PreSmoother.\n");
		return false;
	}
	if(m_postsmoother==NULL)
	{
		UG_LOG("AMG::init_amg(): No PostSmoother selected. Call set_postsmoother(postsmoother) to set a PostSmoother.\n");
		return false;
	}

	// init amghelper for grid printing
	amghelper.positions = &dbg_positions[0];
	amghelper.size = A[0]->num_rows();
	amghelper.parentIndex = parentIndex;
	amghelper.dimension = dbg_dimension;

	UG_LOG("Starting AMG Setup." << endl << endl);

	stopwatch SWwhole;
	SWwhole.start();

	int i=0;
	while(i< max_levels-1)
	{
		SMO[i].setmatrix(A[i]);
		m_presmoothers[i] = m_presmoother->clone();
		m_presmoothers[i]->init(SMO[i]);
		if(m_presmoother == m_postsmoother)
			m_postsmoothers[i] = m_presmoothers[i];
		else
		{
			m_postsmoothers[i] = m_postsmoother->clone();
			m_postsmoothers[i]->init(SMO[i]);
		}

		double L = A[i]->num_rows();

#ifndef LATE_COARSE_SOLVER
		UG_LOG("No Late Coarse Solver\n");
		//if(L < 100 || A[i]->total_num_connections()/(L*L) > 0.5)	break; // abbruch falls density > 50%
		if(L < 1000)	break; // abbruch falls density > 50%
#else
		UG_LOG("Late Coarse Solver\n");
		if(L < 10)	break;
#endif
		//smoother[i].init(*A[i]);

		A[i+1] = new matrix_type;
#ifdef UG_PARALLEL
		A[i+1]->set_storage_type(PST_ADDITIVE);
#endif
		create_AMG_level(*A[i+1], R[i], *A[i], P[i], i);

		i++;
	}

	UG_ASSERT(block_vector_traits< typename vector_type::value_type >::is_static, "dynamic not yet implemented");
	int static_nrUnknowns = block_vector_traits< typename vector_type::value_type >::static_size;

	UG_LOG("Creating level " << i << " (" << A[i]->num_rows() << " nodes, total "
			<< A[i]->num_rows()*static_nrUnknowns << " unknowns)" << endl << "Using Direct Solver on Matrix "
			<< A[i]->num_rows()*static_nrUnknowns << "x" << A[i]->num_rows()*static_nrUnknowns << ". ");

	stopwatch SW; SW.start();
	SMO[i].setmatrix(A[i]);
	m_basesolver->init(SMO[i]);

	UG_LOG("Coarse Solver Setup took " << SW.ms() << "ms." << endl);

	used_levels = i+1;
	UG_LOG("AMG Setup finished. Used Levels: " << used_levels << ". ");
	UG_LOG("AMG Setup took " << SWwhole.ms() << " ms." << endl);

	// calc complexities
	double nnzs=0;
	double totallength=0;
	for(int i=0; i<used_levels; i++)
	{
		nnzs += A[i]->total_num_connections();
		totallength += A[i]->num_rows();
	}

	UG_LOG("Operator Complexity: " << nnzs/A[0]->total_num_connections() << " nodes complexity: "
			<< totallength/A[0]->num_rows() << endl << endl);

	return true;
}


//!
//! amg constructor
template<typename TAlgebra>
amg<TAlgebra>::amg() :
	nu1(2),
	nu2(2),
	gamma(1),

	eps_truncation_of_interpolation(0.3),
	theta(0.3),
	sigma(0.3),

	max_levels(10),
	used_levels(0),
	aggressiveCoarsening(0),
	aggressiveCoarseningNrOfPaths(2),

	m_presmoother(NULL),
	m_postsmoother(NULL),
	m_basesolver(NULL)
{
	vec4 = NULL;
	m_bInited = false;


#ifdef UG_PARALLEL
	UG_ASSERT(pcl::GetNumProcesses() == 1, "AMG currently only for 1 process");
#endif

	//FORCE_CREATION { printCoarsening(0,0); }
}


template<typename TAlgebra>
void amg<TAlgebra>::cleanup()
{
	for(int i=1; i<used_levels-1; i++)
	{
		if(A[i]) delete A[i]; A[i] = NULL;
		if(vec1[i]) delete vec1[i]; vec1[i] = NULL;
		if(vec2[i]) delete vec2[i]; vec2[i] = NULL;
		if(parentIndex[i]) delete [] parentIndex[i]; parentIndex[i] = NULL;
	}

	for(int i=0; i<used_levels; i++)
	{
		if(m_presmoothers[i]) delete m_presmoothers[i]; m_presmoothers[i] = NULL;
		if(m_postsmoothers[i]) delete m_postsmoothers[i]; m_postsmoothers[i] = NULL;
		if(vec3[i]) delete vec3[i]; vec3[i] = NULL;
	}
	if(vec4) { delete[] vec4; vec4 = NULL; }
	A[0] = NULL;
	used_levels = 0;
	m_bInited=false;
}
//!
//! amg destructor
template<typename TAlgebra>
amg<TAlgebra>::~amg()
{
	cleanup();
}

template<typename matrix_type, typename vector_type>
void amg_jacobi_step(const matrix_type &Ah, vector_type &d, vector_type &c, vector_type &corr, number damp)
{
	for(size_t i=0; i<c.size(); i++)
	{
		c[i] = damp*(d[i] / Ah.get_diag(i));
	}
	d -= Ah*c;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// get_correction_and_update_defect:
//------------------------------------

template<typename TAlgebra>
bool amg<TAlgebra>::get_correction_and_update_defect(vector_type &c, vector_type &d, int level)
{
	UG_ASSERT(c.size() == d.size() && c.size() == A[level]->num_rows(),
			"c.size = " << c.size() << ", d.size = " << d.size() << ", A.size = " << A[level]->num_rows() << ": not matching");

	const matrix_type &Ah = *(A[level]);

	if(level == used_levels-1)
	{
		m_basesolver->apply_return_defect(c, d);
		return true;
	}

	vector_type &corr = *vec3[level];

	// presmooth
	// same as setting c.set(0.0).
	m_presmoothers[level]->apply_update_defect(c, d);
	for(int i=1; i < nu1; i++)
	{
		m_presmoothers[level]->apply_update_defect(corr, d);
		c += corr;
	}

	vector_type &cH = *vec1[level+1];
	vector_type &dH = *vec2[level+1];

	cH.set_storage_type(PST_CONSISTENT);

	// restrict defect
	// dH = R[level]*d;
	dH.set_storage_type(PST_ADDITIVE);
	MatMult(dH, 1.0, R[level], d);

	// apply lmgc on coarser nodes
	if(level+1 == used_levels-1)
		get_correction_and_update_defect(cH, dH, level+1);
	else
		for(int i=0; i<gamma; i++)
			get_correction_and_update_defect(cH, dH, level+1);

	// interpolate correction
	// corr = P[level]*cH
	corr.set_storage_type(PST_ADDITIVE);
	MatMult(corr, 1.0, P[level], cH);

	// add coarse grid correction to level correction
	// c += corr;
	corr.change_storage_type(PST_CONSISTENT);
	VecScaleAdd(c, 1.0, c, 1.0, corr);

	//update defect
	// d = d - Ah*corr
	MatMultAdd(d, 1.0, d, -1.0, Ah, corr);

	// postsmooth
	for(int i=0; i < nu2; i++)
	{
		m_postsmoothers[level]->apply_update_defect(corr, d);
		c += corr;
	}

	return true;
}



template<typename TAlgebra>
bool amg<TAlgebra>::get_correction(vector_type &c, const vector_type &const_d)
{

	UG_ASSERT(c.size() == const_d.size() && c.size() == A[0]->num_rows(),
				"c.size = " << c.size() << ", d.size = " << const_d.size() << ", A.size = " << A[0]->num_rows() << ": not matching");

	//int level = 0;
	//const matrix_type &Ah = *(A[level]);


	//////////////////
	{
		if(vec4 == NULL)
		{
			vec4 = new vector_type;
			vec4->create(c.size());

		#ifdef UG_PARALLEL
				// todo: change this for later "right" parallel implementation
				UG_ASSERT(pcl::GetNumProcesses() == 1, "AMG currently only for 1 process");
				vec4->set_storage_type(PST_ADDITIVE);
		#endif
		}

		vector_type &d = *vec4;

		d = const_d;
		return get_correction_and_update_defect(c, d);
	}
	//////////////////

/*
	if(level == used_levels-1)
	{
		m_basesolver.apply_return_defect(c, d);
		return true;
	}

	vector_type &corr = *vec3[level];

	// presmooth
	m_presmoothers[level]->apply(c, const_d);

	if(vec4 == NULL)
	{
		vec4 = new vector_type;
		vec4->create(c.size());

#ifdef UG_PARALLEL
		// todo: change this for later "right" parallel implementation
		UG_ASSERT(pcl::GetNumProcesses() == 1, "AMG currently only for 1 process");
		vec4->set_storage_type(PST_ADDITIVE);
#endif
	}
	else if(vec4->size() != c.size())
		vec4->resize(c.size());

	vector_type &d = *vec4;
#ifdef UG_PARALLEL
	d.set_storage_type(PST_ADDITIVE);
#endif
	// d = const_d - Ah*c
	MatMultAdd(d, -1.0, Ah, c, 1.0, const_d);
	for(int i=1; i < nu1; i++)
	{
		// calc c = B^{-1} b
		m_presmoothers[level]->apply_update_defect(c, d);
	}

	vector_type &cH = *vec1[level+1];
	cH.set_storage_type(PST_CONSISTENT);

	vector_type &dH = *vec2[level+1];

	// restrict defect
	// dH = R[level]*d;
	dH.set_storage_type(PST_ADDITIVE);
	MatMult(dH, 1.0, R[level], d);

	// apply lmgc on coarser nodes
	if(level+1 == used_levels-1)
		get_correction_and_update_defect(cH, dH, level+1);
	else
		for(int i=0; i<gamma; i++)
			get_correction_and_update_defect(cH, dH, level+1);

	// interpolate correction
	// corr = P[level]*cH
	corr.set_storage_type(PST_ADDITIVE);
	MatMult(corr, 1.0, P[level], cH);

	// add coarse grid correction to level correction
	// c += corr;
	corr.change_storage_type(PST_CONSISTENT);
	VecScaleAdd(c, 1.0, c, 1.0, corr);

	//update defect
	// d = d - Ah*corr
	MatMultAdd(d, -1.0, Ah, corr, 1.0, d);

		// postsmooth
	for(int i=0; i < nu2-1; i++)
		m_postsmoothers[level]->apply_update_defect(c, d);

	m_postsmoothers[level]->apply(c, d);
*/
	return true;
}

template<typename TAlgebra>
void amg<TAlgebra>::tostring() const
{
	UG_LOG("AMGPreconditioner.\n");
	UG_LOG("nu1 = " << nu1 << endl);
	UG_LOG("nu2 = " << nu2 << endl);
	UG_LOG("gamma = " << gamma << endl);
	UG_LOG("theta = " << theta << endl);
	UG_LOG("sigma = " << sigma << endl);
	UG_LOG("max levels = " << max_levels << endl);

	if(aggressiveCoarsening)	{UG_LOG("Aggressive Coarsening is on, A" << aggressiveCoarseningNrOfPaths << "-mode." << endl);}
	else						{UG_LOG("no Aggressive Coarsening" << endl);}

	if(m_presmoother) 	{UG_LOG("presmoother is " << m_presmoother->name() << ".\n");}
	else				{UG_LOG("no presmoother set!\n");}

	if(m_postsmoother) 	{UG_LOG("postsmoother is " << m_postsmoother->name() << ".\n");}
	else				{UG_LOG("no postsmoother set!\n");}

	if(m_basesolver)	{UG_LOG("basesolver set\n");}
	else				{UG_LOG("no basesolver set!\n");}
}

} // namespace ug

#endif //  __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_IMPL_H__
