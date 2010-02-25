/*
 *  famg.hpp
 *  flexamg
 *
 *  Created by Arne Nägel, Martin Rupp on 25.02.10.
 *  Copyright 2010 . All rights reserved.
 *
 */
#pragma once

extern int *parentIndex[FAMG_MAX_LEVELS];

#if 0
//#define FAMG_WRITE_MATRICES_PATH "/Users/mrupp/matrices/FAMG_"
#define FAMG_WRITE_MATRICES_MAX (200*200)


//#define AMG_PRINT_GRAPH

#define FAMG_WRITE_COARSENING
#define FAMG_WRITE_GRAPH


#define FAMG_PRINT_COARSENING
#define FAMG_PRINT_P
#define FAMG_PRINT_R
#define FAMG_PRINT_AH

#define FAMG_PRINT_COARSEN_RATINGS
#define FAMG_PRINT_COARSEN

#endif

inline double amg_value(double d)
{
	return d;
}

template<typename T>
inline double amg_value(T &d)
{
	return -d.norm();
}
//
// writeMatrices:
//----------------
//! writes A to pathAndName + "A" + level + ".mat for all levels, also R
//! @param pathAndName	path with name of matrix. for example "/Users/username/matrices/mat"
template<typename Matrix_type, typename Vector_type>
void famg<Matrix_type, Vector_type>::writeMatrices(const char *pathAndName)
{
	// only for small matrices
	if(A[0]->getRows() > 100*100*100*100)
		return;
	
	cout << "writing matrices "; cout.flush();
	string str(pathAndName);
	for(int i=0; i<used_levels-1; i++)
	{
		A[i]->writeToFile((str + "A" + nrstring(i) + ".mat").c_str()); cout << "."; cout.flush();		
		P[i].writeToFile((str + "P" + nrstring(i) + ".mat").c_str()); cout << "."; cout.flush();
		R[i].writeToFile((str + "R" + nrstring(i) + ".mat").c_str()); cout << "."; cout.flush();
	}
	if(used_levels > 0)
		A[used_levels-1]->writeToFile((str + "A" + nrstring(used_levels-1) + ".mat").c_str());
	cout << " finished."; cout.flush();
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateGraph:
//-------------------------
//! Create graph of strong connections from A, calculate ratings in grid[i].rating,
//! build up priority queue PQ. unassigned: nr of nodes to be assigned by coarsening algorihm
//!
//! @param	A			matrix A for which to calculate strong connectivity graph
//! @param graph		the calculated strong connectivity graph of A
//!						graph is afterwards made up of connections from a node i to j if
//!						j has a strong connection to i
//! @param	PQ			maxheap priority queue for sorting of the nodes wrt the rating
//! @param unassigned	nr of nodes which are now to be assigned coarse or fine
template<typename Matrix_type, typename Vector_type> // template<typename conn_matrix> // const conn_matrix &C
void famg<Matrix_type, Vector_type>::CreateGraph(const Matrix_type &A, cgraph &graph, maxheap<nodeinfo> &PQ, int &unassigned)
{

	unassigned = 0;
	for(int i=0; i< A.getLength(); i++)
	{
		graph.init(i);
		if(A[i].isUnconnected())
			continue;
		
		double dmax = 0;
		typename Matrix_type::cRowIterator conn = A.beginRow(i); ++conn;
		for(; !conn.isEnd(); ++conn)
		{
			if((*conn).dValue != 0.0 && amg_value((*conn).dValue) < dmax) 
				dmax = amg_value((*conn).dValue);
		}
		
		conn.rewind(); ++conn; // skip diagonal
		for(; !conn.isEnd(); ++conn)
			if( amg_value((*conn).dValue) < theta * dmax)
				graph.setConnection(i, (*conn).iIndex);				
		
	}
#ifdef FAMG_PRINT_GRAPH
	graph.print();
#endif
	
#ifdef FAMG_WRITE_GRAPH
	graph.writeToFile((string(FAMG_WRITE_MATRICES_PATH) + "G" + nrstring(A.tolevel) + ".mat").c_str(), A.tolevel);
#endif
	
	
	// we need the transpose, since when we set a node coarse, we want
	// all nodes to be fine which can be interpolated by this coarse node
	// graph is afterwards made up of connections from a node i to j if
	// j has a strong connection to i
	graph.transpose();
	for(int i=0; i < A.getLength(); i++)
	{
		if(graph.getNrOfConnections(i) == 0)
			grid[i].setFineDirect();
		else
		{
			//ASSERT2(graph.iNrOfConnections[i] > 0, "node " << i << " has " << graph.iNrOfConnections[i] << " connections?");
			grid[i].rating = graph.getNrOfConnections(i);
			PQ.insertItem(i);
			unassigned++;
		}
	}	
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Coarsen:
//-------------------------
//! Coarsens the graph with ratings of nodes in grid[i].rating, set up in a priority queue PQ
//! @param newIndex		store in newIndex[i] new index of node i in coarse grid (>0, if fine < 0)
//! @param unassigned	nr of nodes to assign
//! @param bIndirect	if true, this is 2nd stage of Aggressive Coarsening, then fine nodes get marker "IndirectFine"
//!						instead of just "fine". Used later in CreateProlongation and CreateIndirectProlongation
//! @param A			matrix A (for debug)
//! @return				returns number of new coarse nodes.
template<typename Matrix_type, typename Vector_type>
int famg<Matrix_type, Vector_type>::Coarsen(cgraph &graph, maxheap<nodeinfo> &PQ, int *newIndex, int unassigned, int &iNrOfCoarse, const Matrix_type &A)
{
#ifdef FAMG_WRITE_COARSENING
	fstream fstr((string(FAMG_WRITE_MATRICES_PATH) + "A" + nrstring(A.tolevel) + ".mat").c_str(), ios::out|ios::app);
	fstr << "c" << endl;
#endif
	// construct coarse grid
	//cout << "construct coarse grid" << endl;	
	// old 749 ms bei 1000
	
	while(unassigned > 0)
	{
		// get Node with best rating
		int best = PQ.removeMax();
		
#ifdef FAMG_PRINT_COARSEN
		cout << endl << "set coarse: " << best << " [" << GetOriginalIndex(A.tolevel, best) << "]. rating " << grid[best].rating  << ". then fine: ";
#endif
#ifdef FAMG_WRITE_COARSENING
		fstr << GetOriginalIndex(A.tolevel, best) << endl;
#endif	
		ASSERT2(!grid[best].isAssigned(), "node " << best << " is already assigned??? (rating = " << grid[best].rating << ", unassigend = " << unassigned << ")");
		
		
		// mark as coarse/assigned
		grid[best].setCoarse();
		newIndex[best] = iNrOfCoarse++;
		
		unassigned--;
		
		
		
		// remove neighbors from PQ, so it wont update		
		for(cgraph::cRowIterator conn (graph, best); !conn.isEnd(); ++conn)
		{
			int indexN = conn();
			//cout << graph.conn[best][i] << " ";			
			if(grid[indexN].isAssigned())
				continue;					
			PQ.remove(indexN);			
		}
		
		// mark neighbors as fine
		//cout << " fine: ";
		
		
		for(cgraph::cRowIterator conn (graph, best); !conn.isEnd(); ++conn)
		{
			int indexN = conn();			
			
#ifdef FAMG_PRINT_COARSEN
			cout << indexN << " [" << GetOriginalIndex(A.tolevel, indexN) << "] "<< " ";
			if(grid[indexN].isAssigned())
				cout << (grid[indexN].isCoarse() ? "(c) " : "(f) ");
#endif
			
			if(grid[indexN].isAssigned())
				continue;
			
			//if(bIndirect) grid[indexN].setFineIndirect();
			//else 
			grid[indexN].setFineDirect();
			
			
			unassigned--;
			
			// increase rating of neighbors of neighbors
			
			for(cgraph::cRowIterator connN (graph, indexN); !connN.isEnd(); ++connN)
			{
				int indexNN = connN();
				// TODO: perhaps we could create a the f-f candidate list here
				
				if(grid[indexNN].isAssigned())
					continue;	
				grid[indexNN].rating++;
				PQ.upheap(indexNN);
			}
		}	
		//coarse.print();
		//cout << "Ranking: " << endl;
		//PQ.print();
	}
	//cout << endl;
	
	ASSERT2(iNrOfCoarse > 0, "no coarse nodes???");
	
	cout << iNrOfCoarse << " Coarse Nodes. ";
	
	// second pass
	//----------------
	cgraph TG;
	TG.createAsTransposeOf(graph);
	
	int nrOfFFCoarseNodes=0;
	vector<bool> marks(graph.getLength());
	
#ifdef FAMG_WRITE_COARSENING
	fstr << "c" << endl;
#endif
	
	for(int i=0; i< graph.getLength(); i++)
	{
		if(grid[i].isCoarse() || A.isUnconnected(i))
			continue;		
		
		// mark coarse nodes interpolating this fine node
		for(cgraph::cRowIterator it(TG, i); !it.isEnd(); ++it)
		{
			if(grid[it()].isCoarse())
				marks[it()] = true;
		}
		
		// prevent strong F-F connections without common Interpolation node
		for(cgraph::cRowIterator it(TG, i); !it.isEnd(); ++it)
		{
			if(grid[it()].isCoarse() || A.isUnconnected(it()))
				continue;			
			
			cgraph::cRowIterator it2(TG, it());
			for(; !it2.isEnd(); ++it2)
			{
				if(grid[it2()].isCoarse() && marks[it2()])
					break;
			}	
			
			if(it2.isEnd())
			{
				// TODO: calculate ratings, add to PQ and do coarsening again on those candidates
				// rating of a fine node = nr of f-f pairs which this node is adjacent to BOTH
				// that is: all common fine node neighbors of i and it() get rating++.
				// problem: updating
#ifdef FAMG_WRITE_COARSENING
				fstr << GetOriginalIndex(A.tolevel, i) << endl; // << endl << GetOriginalIndex(A.tolevel, it()) << endl;
#endif
				//cout << endl << "prevent F-F-connection between (2) " << i << "[" << GetOriginalIndex(A.tolevel, i) << "] and " << it() << "[" << GetOriginalIndex(A.tolevel, it()) << "], setting " << i << " coarse.";				
				grid[i].setCoarse();
				newIndex[i] = iNrOfCoarse++;
				
				nrOfFFCoarseNodes++;
				break;
			}
		}
		
		// remove marks
		for(cgraph::cRowIterator it(graph, i); !it.isEnd(); ++it)
		{
			if(grid[it()].isCoarse())
				marks[it()] = false;
		}		
	}	
	
	if(nrOfFFCoarseNodes)
		cout << "F-F prevention, now " << iNrOfCoarse << " coarse nodes." << endl;
	return iNrOfCoarse;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateProlongation:
//-------------------------
//! Calculates Prolongation P with Matrix_type A and coarse/fine markers in grid[i].isFine/isCoarse by direct interpolation
//! uses nodeinfo *grid.
//!	@param P				Matrix P: here goes the calculated prolongation
//! @param	A				Matrix A: matrix for which to calculate prolongation on next level
//! @param	newIndex		newIndex of coarse Node i on next coarser level
//! @param	iNrOfCoarse		nr of coarse nodes on this level
template<typename Matrix_type, typename Vector_type>
void famg<Matrix_type, Vector_type>::CreateProlongation(SparseMatrix<double> &P, const Matrix_type &A, int *newIndex, int iNrOfCoarse, int &unassigned)
{	
	P.create(A.getLength(), iNrOfCoarse);
	P.fromlevel = A.fromlevel+1;
	P.tolevel = A.tolevel;
	P.name = "AMG:P";
	vector<SparseMatrix<double>::connection> con(255);
	SparseMatrix<double>::connection c;
	// DIRECT INTERPOLATION
	unassigned=0;
	
	for(int i=0; i < A.getLength(); i++)
	{
		if(grid[i].isCoarse())
		{
			// a coarse node
			//P[i].initWithoutDiag();
			SparseMatrix<double>::connection con; 
			con.iIndex = newIndex[i];  ASSERT1(newIndex[i] != -1);
			con.dValue = 1.0;
			P.setMatrixRow(i, &con, 1);
		}
		else if(A[i].isUnconnected())
		{
			//P[i].initWithoutDiag(); // boundary values need not to be prolongated
		}
		else if(grid[i].isFineDirect())
		{	
			// a non-interpolated fine node. calculate interpolation weights
			// see amg::CreateProlongation in amg.hpp
			
			
			vector<int> indices;
			blockDenseMatrix<double, variableStorage> SM(indices.size(), indices.size());
			
			// get Neighborhood of depth 2
			A.getNeighborhood(i, 2, indices);
			
			SM.setSize(indices.size(), indices.size());

			// get corresponding submatrix
			A.get(SM, indices, indices);

			// find index of i in local indices.
			int ind = find(indices.begin(), indices.end(), i) - indices.begin();
			
			// traverse off-diag entries of A			
			typename Matrix_type::cRowIterator conn = A.beginRow(i); 
			++conn; // skip diag			
			for(; !conn.isEnd(); ++conn)
			{
				// elements via (*conn).dValue and (*conn).iIndex.
				
			}
			
			
			// traverse graph
			for(cgraph::cRowIterator it(graph, i); !it.isEnd(); ++it)
			{
				int con = it();
			}
			
			
			// connections hinzufügen
			if(con.size() > 0)
			{				
				P.setMatrixRow(i, &con[0], con.size());	
			}
			else
			{
				unassigned++;
				grid[i].setFineIndirect();
			}
		}
		else
		{
			unassigned++;
			ASSERT2(aggressiveCoarsening != 0, "no aggressive Coarsening but node " << i << " is fine and indirect??");
		}		
	}	
	
	if(unassigned)
		cout << "Pass 1: " << unassigned << " left. ";
}

}	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// createAMGLevel:
//-------------------------
//! create AMG matrix R, P, and AH = R A P
//! @param AH
//! @param R
//! @param A
//! @param P
//! @param level
template<typename Matrix_type, typename Vector_type>
void famg<Matrix_type, Vector_type>::createAMGLevel(Matrix_type &AH, SparseMatrix<double> &R, const Matrix_type &A, SparseMatrix<double> &P, int level)
{
	bool bTiming=true;
	cout << "Creating level " << level << ". (" << A.getLength() << " nodes)" << endl;
	stopwatch SW;
	stopwatch SWwhole; SWwhole.start();
	
	// nodeinfo: infos zu den einzelnen knoten von A, verwaltung der ratings etc
	// entries: speicherung der strong connections, nur wohin
	// entriesvalues: speicherung der WERTE der connections.
	// (Werte werden bei Konstruktion des Grobgitters nicht benötigt)
	
	// todo: check for isolated condition
	
	maxheap<nodeinfo> PQ(A.getLength(), grid);
	
	int unassigned = A.getLength();
	
	int *newIndex = new int[A.getLength()];
	memset(newIndex, -1, sizeof(int)*A.getLength());
	
	int *posInConnections = new int[A.getLength()];
	memset(posInConnections, -1, sizeof(int)*A.getLength());
	
	cgraph graph(A.getLength(), A.getTotalNrOfConnections()*10);
	
	
	// build graph
	/////////////////////////////////////////
	
	cout << "building graph... "; cout.flush();
	if(bTiming) SW.start();
	// CreateGraph(C, graph, PQ, unassigned);
	CreateGraph(A, graph, PQ, unassigned);	
	if(bTiming) SW.printTimeDiff();
	
#ifdef FAMG_PRINT_COARSEN_RATINGS
	for(int i=0; i<A.getLength(); i++)
	{
		cout << i << " (" << GetPosForIndexAtLevel(i, level).x << " " << GetPosForIndexAtLevel(i, level).y << ") ";
		grid[i].print();
	}
#endif
	//PQ.print();
	
	// Coarsen
	/////////////////////////////////////////
	cout << "coarsening... "; cout.flush();	
	
	if(bTiming) SW.start();
	int iNrOfCoarse = 0;
	Coarsen(graph, PQ, newIndex, unassigned, iNrOfCoarse, A);	
	
	if(bTiming) { SW.printTimeDiff();}
	if(bTiming) cout << endl;
	
	// create vectors for AMG multigrid
	/////////////////////////////////////////
	
	vec1[level+1] = new Vector_type (iNrOfCoarse, "AMG:tempvec 1");
	vec1[level+1]->level = level+1;
	vec2[level+1] = new Vector_type (iNrOfCoarse, "AMG:tempvec 2");
	vec2[level+1]->level = level+1;
	cout << "created vec1 on level" << level +1 << endl;
	
	// set size for variable sized blockvectors
	for(int i=0; i<A.getLength(); i++)
		if(grid[i].isCoarse())
		{			
			setSize((*vec1[level+1])[newIndex[i]], getRows(A[i][0].dValue));
			setSize((*vec2[level+1])[newIndex[i]], getRows(A[i][0].dValue));
		}
	
	// set parentindex for debugging
	parentIndex[level+1] = new int[iNrOfCoarse];
	for(int i=0; i<A.getLength(); i++)
		if(grid[i].isCoarse())
			parentIndex[level+1][ newIndex[i] ] = i;
	
	/*for(int i=0; i<A.getLength(); i++)
	 {
	 cout << i << " (" << GetPosForIndexAtLevel(i, level).x << " " << GetPosForIndexAtLevel(i, level).y << ") ";
	 grid[i].print();
	 }//*/
	
#ifdef FAMG_PRINT_COARSENING
	printCoarsening(level);
#endif	
	
	// construct prolongation P = I_{2h->h}
	/////////////////////////////////////////
	
	cout << "prolongation... "; cout.flush();
	
	unassigned = 0;
	
	P.fromlevel = level+1;
	P.tolevel = level;
	
	if(bTiming) SW.start();
	CreateProlongation(P, A, newIndex, iNrOfCoarse, unassigned);
	if(unassigned > 0)
		CreateIndirectProlongation(P, A, newIndex, posInConnections, unassigned);
	
	if(bTiming) SW.printTimeDiff();
	
#ifdef FAMG_PRINT_P
	P.print();
#endif
	
	// construct prolongation R = I_{h->2h}
	/////////////////////////////////////////
	
	cout << "restriction... "; cout.flush();
	if(bTiming) SW.start();
	// construct restriction R = I_{h -> 2h}		
	R.createAsTransposeOf(P); // already finished
	R.name = "AMG:R";
	//R.print("R");
	if(bTiming) SW.printTimeDiff();	
	
#ifdef FAMG_PRINT_R
	R.print();
#endif
	
	// create Galerkin product
	/////////////////////////////////////////	
	
	cout << "galerkin product... "; cout.flush();
	if(bTiming) SW.start();
	
	// AH = R A P
	AH.createAsMultiplyOf(R, A, P, posInConnections);
	
	if(bTiming) SW.printTimeDiff();
	AH.name = "AMG:A";
	AH.fromlevel = level+1;
	AH.tolevel = level+1;
	AH.finish();
	
#ifdef FAMG_PRINT_AH
	AH.print();
#endif
	
	// finish
	/////////////////////////////////////////	
	
	//AH.print("AH");
	int nnz = AH.getTotalNrOfConnections();
	cout << "AH: nnz: " << nnz << " Density: " << double(nnz)/(double(AH.getLength())*double(AH.getLength()))*100.0 << "% nnz/n: " << nnz/(double)AH.getLength() << endl;
	cout << "Coarsening rate: " << (100.0*AH.getLength())/(A.getLength()) << "%" << endl;
	
	cout << " level "; SWwhole.printTimeDiff();  cout << endl;
	cout << endl;
	cout.flush();
	
	cout << A << endl;
	cout << P << endl;
	cout << R << endl;
	
	
#ifdef FAMG_WRITE_MATRICES_PATH
	if(this->A[0]->getLength() < FAMG_WRITE_MATRICES_MAX)
	{
		cout << "write matrices";
		P.writeToFile((string(FAMG_WRITE_MATRICES_PATH) + "P" + nrstring(level) + ".mat").c_str()); cout << "."; cout.flush();
		R.writeToFile((string(FAMG_WRITE_MATRICES_PATH) + "R" + nrstring(level) + ".mat").c_str()); cout << "."; cout.flush();
		AH.writeToFile((string(FAMG_WRITE_MATRICES_PATH) + "A" + nrstring(level+1) + ".mat").c_str()); cout << "."; cout.flush();
		cout << " done." << endl;
	}
#endif
	
	//P.print();
	
	cout << AH << endl;
	cout << endl;//*/
	
	delete[] posInConnections;
	delete[] newIndex;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// famg<Matrix_type, Vector_type>::onlyOneLevel
//----------------
//! for testing. creates only one AMG level without direct solvers
//! @param A	matrix A.
template<typename Matrix_type, typename Vector_type>
bool famg<Matrix_type, Vector_type>::onlyOneLevel(const Matrix_type& A_)
{
	used_levels = 2;
	const Matrix_type *pA = &A_;
	grid = new nodeinfo[pA->getLength()];
	A[0] = const_cast<Matrix_type*> (pA);
	
#ifdef FAMG_WRITE_MATRICES_PATH	
	if(A[0]->getLength() < FAMG_WRITE_MATRICES_MAX)
	{
		cout << "write matrix A...";
		A[0]->writeToFile((string(FAMG_WRITE_MATRICES_PATH) + "A" + nrstring(0) + ".mat").c_str());
		cout << "done." << endl; cout.flush();
	}
#endif
	
	
	int i=0;
	A[i+1] = new Matrix_type();
	createAMGLevel(*A[i+1], R[i], *A[i], P[i], i);
	
	//	A[1]->print();
	return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// famg<Matrix_type, Vector_type>::init
//----------------
//! creates MG Hierachy for with Matrix_type A and temporary vectors for higher levels
//! @param A	matrix A.
template<typename Matrix_type, typename Vector_type>
bool famg<Matrix_type, Vector_type>::init(const Matrix_type& A_)
{
	cout << "Starting AMG Setup." << endl << endl;
	stopwatch SW;
	const Matrix_type *pA = &A_;
	grid = new nodeinfo[pA->getLength()];
	A[0] = const_cast<Matrix_type*> (pA);
	
	A[0]->fromlevel = 0;
	A[0]->tolevel = 0;
	
	
#ifdef FAMG_WRITE_MATRICES_PATH	
	if(A[0]->getLength() < FAMG_WRITE_MATRICES_MAX)
	{
		cout << "write matrix A...";
		A[0]->writeToFile((string(FAMG_WRITE_MATRICES_PATH) + "A" + nrstring(0) + ".mat").c_str());
		cout << "done." << endl; cout.flush();
	}
#endif
	
	int i=0;
	while(i< max_levels-1)
	{
		
		double L = A[i]->getLength();
#if 1
		if(L < 100 || A[i]->getTotalNrOfConnections()/(L*L) > 0.5)	break; // abbruch falls density > 50%
#else
		if(L < 10)	break;
#endif
		smoother[i].init(*A[i]);
		
		A[i+1] = new Matrix_type();
		createAMGLevel(*A[i+1], R[i], *A[i], P[i], i);
		
		vec3[i] = new Vector_type (A[i]->getLength(), "AMG:tempvec 3");
		vec3[i]->level = i;		
		i++;
	}
	
	int nrOfUnknowns = vec_traits< typename Vector_type::entry_type >::nrOfUnknowns;
	cout << "Creating level " << i << " (" << A[i]->getLength() << " nodes, total " << A[i]->getLength()*nrOfUnknowns << " unknowns)" << endl << "Using Direct Solver on Matrix " 
	<< A[i]->getLength()*nrOfUnknowns << "x" << A[i]->getLength()*nrOfUnknowns << ". ";
	stopwatch SW2; SW2.start();
	coarseSolver.create(*A[i]);	
	SW2.printTimeDiff();
	cout << endl << endl;
	
	used_levels = i+1;
	cout << "AMG Setup finished. Used Levels: " << used_levels << ". ";
	SW.printTimeDiff();
	
	// calc complexities
	double nnzs=0;
	double totallength=0;
	for(int i=0; i<used_levels; i++)
	{
		nnzs += A[i]->getTotalNrOfConnections();
		totallength += A[i]->getLength();
	}	
	
	cout << "Operator Complexity: " << nnzs/A[0]->getTotalNrOfConnections() << " grid complexity: " << totallength/A[0]->getLength() << endl << endl;
	
	// cleanup
	delete [] grid;
	
	return true;
}

//!
//! amg constructor
template<typename Matrix_type, typename Vector_type>
famg<Matrix_type, Vector_type>::famg()
{
	used_levels = 0;
	max_levels = 10;
	
	nu1 = 2;
	nu2 = 2;
	gamma = 1;
	
	eps_truncation_of_interpolation = 0.3; // no truncation (or 0.2).
	
	sigma = 0.3;
	theta = 0.5;
	
	//	if(never_happens) printCoarsening(0,0);
}

//!
//! amg destructor
template<typename Matrix_type, typename Vector_type>
famg<Matrix_type, Vector_type>::~amg()
{
	for(int i=1; i<used_levels-1; i++)
	{
		delete A[i];
		delete vec1[i];
		delete vec2[i];
		delete vec3[i-1];
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MGCycle:
//-------------------------
//! calculates one AMG cycle
//! @param		x		
//! @param	b	
//! @param	level	
template<typename Matrix_type, typename Vector_type>
double famg<Matrix_type, Vector_type>::MGCycle(Vector_type &x, const Vector_type &b, int level)
{	
	ASSERT2(x.getLength() == b.getLength() && b.getLength() == A[level]->getLength(), 
			"x (" << x << "), b (" << b << "), or A (" << *A[level] << ") have different length");
	
	const Matrix_type &Ah = *(A[level]);
	
	if(level == used_levels-1)
	{
		coarseSolver.solve(b, x);
		return 0.1e-14;
	}
	
	for(int i=0; i < nu1; i++)
		smoother[level].iterate(x, b);
	Vector_type &r = *vec3[level];
	
	r = b - Ah*x;
	
	Vector_type &rH = *vec1[level+1];
	Vector_type &eH = *vec2[level+1];
	
	rH = R[level]*r;
	
	eH = 0.0;
	
	if(level+1 == used_levels-1)
		MGCycle(eH, rH, level+1);
	else
		for(int i=0; i<gamma; i++)
			MGCycle(eH, rH, level+1);
	
	x += P[level]*eH;
	
	double res;
	//if(level != 0)
	for(int i=0; i < nu2; i++)
		res = smoother[level].iterate(x, b);
	
	//cout << level << " postsmooth " << norm(b - Ah*x) << endl;
	
	return res;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// amgTestLevel:
//-------------------------
//! tests one AMG level by smoothing x, then doing MGCycle down, smoothing, then calc reduction
//! with this function you can tests where your MG hangs
//! @param	x	a 
//! @param	b	
//! @param	level	
template<typename M
template<typename Matrix_type, typename Vector_type>
void famg<Matrix_type, Vector_type>::amgTestLevel(Vector_type &x, const Vector_type &b, int level)
{
	cout.flush(); 
	const Matrix_type &Ah = *(A[level]);
	
	if(level == used_levels-1)
	{
		cout << endl;
		cout << "[" << level << "]." << endl;
		cout << "Coarse Solver. "<< endl;
		coarseSolver.solve(b, x);
		cout << "res: " << norm(b-Ah*x) << endl;
		return;
	}
	
	double pre1 = norm(b - Ah*x);
	
	for(int i=0; i < nu1; i++)
		smoother[level].iterate(x, b);
	Vector_type &r = *vec3[level];
	
	r = b - Ah*x;
	double presmoothreduction = norm(r)/pre1;
	//cout << setw(2*level) << "presmooth reduction: " << norm(r)/pre1 << endl;
	
	Vector_type &rH = *vec1[level+1];
	Vector_type &eH = *vec2[level+1];
	
	rH = R[level]*r; 
	
	eH = 0.0;
	
	if(level+1 == used_levels-1)
		MGCycle(eH, rH, level+1);
	else
		//for(int i=0; i<2; i++)
		MGCycle(eH, rH, level+1);
	
	x += P[level]*eH; //interpolate(eH, level);	
	double pre2 = norm(b - Ah*x);
	
	double res;
	for(int i=0; i < nu2; i++)
		res = smoother[level].iterate(x, b);
	
	double post = norm(b - Ah*x);
	
	//cout << setw(2*level) << "postsmooth reduction: " << post/pre2 << endl;
	//cout << setw(2*level) << "level reduction: " << post/pre1 << endl;
	
	cout << endl;
	cout << "[" << level << "]" << endl;
	cout << "prered: " << presmoothreduction << " postred: " << post/pre2
	<< " levelred: " << post/pre1;
	
	if(post < 1e-10 || post/pre1 < 0.3)
		cout << " -> level seems to be OK." << endl;
	else
		cout << " -> LEVEL BROKEN!" << endl;	
	
	cout << "norm(b-Ah*x) = " << post << endl;	
	
	vector<sortStruct<double> > inner;
	vector<sortStruct<double> > boundary;
	
	for(int i=0; i<r.getLength(); i++)
	{
		sortStruct<double> s;
		s.sortValue = -mnorm(r[i]);
		s.index = i;
		if(level == 0 && Ah.isCloseToBoundary(i, 2))
			boundary.push_back(s);
		else
			inner.push_back(s);
	}
	sort(inner.begin(), inner.end());
	sort(boundary.begin(), boundary.end());
	
	if(boundary.size() > 0 && -boundary[0].sortValue > 1e-7)
	{
		cout << "big error boundary: " << endl;
		for(int i=0; i<10 && i<boundary.size(); i++)
			cout << boundary[i].index << "[" << GetOriginalIndex(level, boundary[i].index) << "]: " << -boundary[i].sortValue << endl;
	}
	if(inner.size() > 0 && -inner[0].sortValue > 1e-7)
	{
		cout << "big error inside: " << endl;
		for(int i=0; i<10 && i<inner.size(); i++)
			cout << inner[i].index << "[" << GetOriginalIndex(level, inner[i].index) << "]: " << -inner[i].sortValue << endl;
		
	}
	
	amgTestLevel(eH, rH, level+1);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// amgTest:
//-------------------------
//! tests all AMG level calling amgTestLevel for all levels.
//! Does some iterations of AMG MG before to ensure that we observe convergence rates at the end
template<typename Matrix_type, typename Vector_type>
void famg<Matrix_type, Vector_type>::amgTest(const Matrix_type& A_, Vector_type &vx, const Vector_type &b)
{
	init(A_);
	
	cout << "++++++ ANALYZING ++++++"<< endl;
	cout << "+++++++++++++++++++++++" << endl;
	
	int i;
	for(i=0; i<10 && (absmax(b - A_*vx) > 1e-6); i++)
		iterate(vx, b);
	cout << endl << i << " presteps. norm: " << norm(b-A_*vx) << endl;
	
	amgTestLevel(vx, b, 0);		
}
