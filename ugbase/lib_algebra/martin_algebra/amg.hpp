/*
 *  amg.hpp
 *  flexamg
 *
 *  Created by Martin Rupp on 03.12.09.
 *  Copyright 2009 . All rights reserved.
 *
 */
#pragma once
//static const double theta = 0.3;

extern int *parentIndex[AMG_MAX_LEVELS];

#define AMG_WRITE_MATRICES_PATH "/Users/mrupp/matrices/AMG_"
#define AMG_WRITE_MATRICES_MAX 10003

#define AMG_PRINT_COARSEN


#if 0
#define AMG_PRINT_COARSENING

#define AMG_PRINT_AH
#define AMG_PRINT_P
#define AMG_PRINT_R

#define AMG_PRINT_GRAPH

#endif

//
// writeMatrices:
//----------------
//! writes A to pathAndName + "A" + level + ".mat for all levels, also R
//! @param pathAndName	path with name of matrix. for example "/Users/username/matrices/mat"
template<typename Matrix_type, typename Vector_type>
void amg<Matrix_type, Vector_type>::writeMatrices(const char *pathAndName)
{
	// only for small matrices
	if(A[0]->getRows() > 1000*1000)
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



// printCoarsening:
//----------------
//! Debug output. Writes position of Coarse nodes in coarse<level>.dat, and fine in fine<level>.dat for display in gnuplot
//! @param level	level which is to be printed
template<typename Matrix_type, typename Vector_type>
void amg<Matrix_type, Vector_type>::printCoarsening(int level)
{  
	fstream fcoarse((string("/Users/mrupp/matrices/coarse") + nrstring(level) + ".dat").c_str(), ios::out);	
	fstream ffine  ((string("/Users/mrupp/matrices/fine") + nrstring(level) + ".dat").c_str(), ios::out);	
	int n = A[level]->getRows();
	
	for(int i=0; i < n; i++)
	{
		postype pos = GetPosForIndexAtLevel(i, level);
		if(grid[i].isCoarse())
			fcoarse << pos.x << " " << pos.y << " " << endl;
		else
			ffine << pos.x << " " << pos.y << " " << endl;
	}
  	
	/////////////
	
	fstream file((string("/Users/mrupp/matrices/coarsening") + nrstring(level) + ".mat").c_str(), ios::out);	
	writePosToStream(file);
	file << 0 << endl;
	for(int i=0; i < n; i++)
	{
		if(!grid[i].isCoarse())
		{
			int org =  GetOriginalIndex(level, i);
			file << org << " " << org << " " << 1.0 << endl;
		}
	}	
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// createGalerkinMatrix:
//-------------------------
// Calculates AH = R A P. posInConnections only needed for speedup (has to be -1 forall i).
template<typename Matrix_type, typename Vector_type>
void amg<Matrix_type, Vector_type>::createGalerkinMatrix(Matrix_type &AH, const SparseMatrix<double> &R, const Matrix_type &A, const SparseMatrix<double> &P, int *posInConnections)
{
	cout << endl << "amg::createGalerkinMatrix: DEPRECATED FUNCTION" << endl;
	int n = R.getLength();
	
	// speedup with array posInConnections, needs n memory (also used in CreateGraph2 and CreateIndirectProlongation).
	// 1000x1000 ninepoint: old version: 2444 ms, new 1200 ms, 1000 ms, 900 ms, 950 ms
	
	// posInConnections[i]: index in the connections for current row.
	// has to be -1 for all nodes	
	
	vector<typename Matrix_type::connection > con(255);
	double r, p;
	entry_type ra;
	
	AH.create(n, n);	
	AH.fromlevel = P.fromlevel;
	AH.tolevel = R.tolevel;
	for(int i=0; i < n; i++)
	{
		//cout << endl << "Node " << i << " is connected to ";
		R.prefetch(i+2);
		
		// we want to have the diagonal first:
		posInConnections[i] = 0;
		con.clear();
		typename Matrix_type::connection c;
		c.iIndex = i;
		c.dValue = 0.0;
		con.push_back(c);
		
		for(typename SparseMatrix<double>::cRowIterator itR(R, i); !itR.isEnd(); ++itR)
		{
			r = (*itR).dValue;
			if(r == 0.0) continue;
			
			for(typename Matrix_type::cRowIterator itA(A, (*itR).iIndex); !itA.isEnd(); ++itA)
			{
				ra = (*itA).dValue * r;
				if(ra == 0.0) continue;
				for(typename SparseMatrix<double>::cRowIterator itP(P, (*itA).iIndex);  !itP.isEnd(); ++itP)
				{
					p = (*itP).dValue;
					if(p == 0.0) continue;
					int indexTo = (*itP).iIndex;					
					
					if(posInConnections[indexTo] == -1)
					{
						// we havent visited node <indexTo>
						// so we need to add a Connection to the row
						// save the index of the connection in the row
						posInConnections[indexTo] = con.size();
						c.iIndex = indexTo;
						c.dValue = ra*p;
						con.push_back(c);				
					}
					else
					{
						// we have visited this node before,
						// so we know the index of the connection
						// -> add r*a*p
						//TODO 
						con[posInConnections[indexTo]].dValue += ra*p;
					}
					
				}
			}
		}
		// set matrix_type Row in AH
		AH.setMatrixRow(i, &con[0], con.size());
		
		// reset posInConnections to -1
		for(int j=0; j<con.size(); j++) posInConnections[con[j].iIndex] = -1;
	}
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
void amg<Matrix_type, Vector_type>::CreateGraph(const Matrix_type &A, cgraph &graph, maxheap<nodeinfo> &PQ, int &unassigned)
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
			if((*conn).dValue != 0.0 && dmax < mnorm((*conn).dValue)) 
				dmax = mnorm((*conn).dValue);
		}
		
		
		grid[i].rating = 0;
		conn.rewind(); ++conn; // skip diagonal
		for(; !conn.isEnd(); ++conn)
			if( mnorm((*conn).dValue) >= theta * dmax)
				graph.setConnection(i, (*conn).iIndex);				
		
	}
#ifdef AMG_PRINT_GRAPH
	graph.print();
#endif
	
	// we need the transpose, since when we set a node coarse, we want
	// all nodes to be fine which can be interpolated by this coarse node
	// graph is afterwards made up of connections from a node i to j if
	// j has a strong connection to i
	graph.transpose();
	for(int i=0; i < A.getLength(); i++)
	{
		if(A[i].isUnconnected())
			grid[i].setFine();
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
// CreateGraph2:
//-------------------------
//! Create the graph2, which consists of coarse nodes in graph.
//! two coarse nodes a and b in graph2 are connected, if
//! - there exist at least 2 ways of length 2 from a to b (A2-Coarsening)
//! - there exist a ways of length 2 from a to b (A1-Coarsening)
//! because of the coarsening process, every coarse node has only fine neighbors in graph,
//! that means those ways are from a coarse node over a fine node to a coarse node.
//!
//! @param	graph			old graph of normal strong connectivity (from CreateGraph)
//! @param graph2			new graph of "distance-2-strong-connectivity"
//! @param PQ				maxheap priority queue for sorting of the nodes wrt the rating
//! @param unassigned		nr of nodes which are now to be assigned coarse or fine
//! @param posInConnections		array of size graph.size for speedup of neighbor-neighbor-calculation inited with -1.
//! @note could be faster with using std::map instead of vector<int> and posInConnections
template<typename Matrix_type, typename Vector_type>
void amg<Matrix_type, Vector_type>::CreateGraph2(cgraph &graph, cgraph &graph2, maxheap<nodeinfo> &PQ, int &unassigned, int *posInConnections)
{		
	vector<int> connection(255);
	vector<int> nrOfPaths(255);
	
	PQ.reset();
	//graph.print();
	unassigned=0;
	for(int b=0; b < graph.getLength(); b++)
	{
		graph2.init(b);
		if(grid[b].isFine())
			continue;
		
		connection.clear(); 
		nrOfPaths.clear();
		// first calculate all nodes reachable with paths of length 2
		
		// ! i is coarse -> has only fine neighbors
		for(cgraph::cRowIterator conn (graph, b); !conn.isEnd(); ++conn)
		{
			int indexN = conn();					
			for(cgraph::cRowIterator connN (graph, indexN); !connN.isEnd(); ++connN)
			{
				int indexNN = connN();
				
				if(indexNN == b || grid[indexNN].isFine())
					continue;
				int pos = posInConnections[indexNN];
				if(pos == -1)
				{
					// never reached node indexNN from b, init.
					pos = posInConnections[indexNN]= connection.size();					
					connection.push_back(indexNN);
					nrOfPaths.push_back(1);					
				}
				else					
					nrOfPaths[pos]++;
			}
		}
		
		// then sort out those which were reached #aggressiveCoarseningNrOfPaths (2 or 1) times
		grid[b].rating = 0;
		for(int i=0; i<connection.size(); i++)
		{
			if(nrOfPaths[i] >= aggressiveCoarseningNrOfPaths)
			{
				// add connection b -> node
				graph2.setConnection(b, connection[i]);
				// increase rating of b
				grid[b].rating ++;
			}
			// reset posInConnections for further use
			posInConnections[connection[i]] = -1;
		}
		
		// add node with rating > 0 to priority queue
		if(grid[b].rating > 0)
		{
			PQ.insertItem(b);
			unassigned++;
		}
		else
			grid[b].setCoarse();
		
	}	
	
	//cout << endl << endl;
	//graph2.print();
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
int amg<Matrix_type, Vector_type>::Coarsen(cgraph &graph, maxheap<nodeinfo> &PQ, int *newIndex, int unassigned, bool bIndirect, const Matrix_type &A)
{
	int iNrOfCoarse	=0;
	// construct coarse grid
	//cout << "construct coarse grid" << endl;	
	// old 749 ms bei 1000
	
	while(unassigned > 0)
	{
		// get Node with best rating
		int best = PQ.removeMax();
		
#ifdef AMG_PRINT_COARSEN
		cout << endl << "set coarse: " << best << " [" << GetOriginalIndex(A.tolevel, best) << "]. rating " << grid[best].rating  << ". then fine: ";
#endif
		
		ASSERT2(!grid[best].isAssigned(), "node " << best << " is already assigned??? (rating = " << grid[best].rating << ", unassigend = " << unassigned << ")");
		
		newIndex[best] = iNrOfCoarse++;
		
		// mark as coarse/assigned
		grid[best].setCoarse();
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
			
#ifdef AMG_PRINT_COARSEN
			cout << indexN << " [" << GetOriginalIndex(A.tolevel, indexN) << "] "<< " ";
			if(grid[indexN].isAssigned())
				cout << grid[indexN].isCoarse() ? "(c) " : "(f) ";
#endif
			
			if(grid[indexN].isAssigned())
				continue;
			
			if(bIndirect) grid[indexN].setFineIndirect();
			else grid[indexN].setFine();
			
			
			unassigned--;
			
			// increase rating of neighbors of neighbors
			
			for(cgraph::cRowIterator connN (graph, indexN); !connN.isEnd(); ++connN)
			{
				int indexNN = connN();
				/*if(grid[indexNN].isFine())
				 {
				 // prevent F-F connections without common Interpolation node
				 int nrOfConnectionsNN = graph.iNrOfConnections[indexNN];
				 
				 int k;
				 for(k=0; k<nrOfConnectionsNN; k++)
				 {
				 int indexNNN = graph.conn[indexNN][k];						
				 if(grid[indexNNN].isCoarse() && indexNNN == best)
				 break;
				 }					
				 if(k == nrOfConnectionsNN)						
				 {
				 cout << "prevent F-F connection between " << indexN << " and " << indexNN << endl;
				 grid[indexN].setCoarse();
				 newIndex[indexN] = iNrOfCoarse++;
				 continue;
				 }
				 }*/
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
	
	//return iNrOfCoarse;
	
	// second pass
	//----------------
	// seems to work but doesnt help with convergence rates on complicated geometries???	
	cgraph TG;
	TG.createAsTransposeOf(graph);
	
	vector<bool> marks(graph.getLength());
	for(int i=0; i< graph.getLength(); i++)
	{
		if(!grid[i].isFine())
			continue;
		
		
		// mark coarse nodes interpolating this fine node
		for(cgraph::cRowIterator it(graph, i); !it.isEnd(); ++it)
		{
			if(grid[it()].isCoarse())
				marks[it()] = true;
		}
		
		
		// prevent strong F-F connections without common Interpolation node
		for(cgraph::cRowIterator it(graph, i); !it.isEnd(); ++it)
		{
			if(!grid[it()].isFine())
				continue;
			
			
			cgraph::cRowIterator it2(TG, it());
			for(; !it2.isEnd(); ++it2)
			{
				if(grid[it2()].isCoarse() && marks[it2()])
					break;
			}	
			
			if(it2.isEnd())
			{
				//cout << "prevent F-F-connection between " << i << " and " << it() << endl;				
				grid[i].setCoarse();
				newIndex[i] = iNrOfCoarse++;
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
void amg<Matrix_type, Vector_type>::CreateProlongation(SparseMatrix<double> &P, const Matrix_type &A, int *newIndex, int iNrOfCoarse)
{	
	P.create(A.getLength(), iNrOfCoarse);
	P.fromlevel = A.fromlevel+1;
	P.tolevel = A.tolevel;
	P.name = "AMG:P";
	vector<SparseMatrix<double>::connection> con(255);
	SparseMatrix<double>::connection c;
	// DIRECT INTERPOLATION
	
	for(int i=0; i < A.getLength(); i++)
	{
		if(grid[i].isCoarse())
		{
			// a coarse node
			//P[i].initWithoutDiag();
			SparseMatrix<double>::connection con; 
			con.iIndex = newIndex[i]; 
			con.dValue = 1.0;
			P.setMatrixRow(i, &con, 1);
		}
		else if(A[i].isUnconnected())
		{
			//P[i].initWithoutDiag(); // boundary values need not to be prolongated
		}
		else if(!grid[i].isFineIndirect())
		{	
			// a non-interpolated fine node. calculate interpolation weights
			
			const double eps_truncate = 0.2;	// Ruge/Stuebe A.7.2.4 truncation of interpolation
			
			// calc min off-diag-entry
			double dmax = 0, connValue, maxConnValue = 0;
			double sumNeighbors =0, sumInterpolatory=0;
			typename Matrix_type::cRowIterator conn = A.beginRow(i); ++conn; // skip diag
			
			for(; !conn.isEnd(); ++conn)
			{
				connValue = mnorm((*conn).dValue);
				sumNeighbors += -connValue;
				
				if(dmax < connValue) 
					dmax = connValue;
				if(!grid[(*conn).iIndex].isFine() && maxConnValue < connValue)
					maxConnValue = connValue;
				
			}
			
			double barrier = min(theta*dmax, eps_truncate*maxConnValue);
			con.clear();			
			
			conn.rewind(); ++conn; // skip diagonal
			for(; !conn.isEnd(); ++conn)
			{
				connValue = -mnorm((*conn).dValue);								
				if(-connValue < barrier || grid[(*conn).iIndex].isFine())
					continue;
				c.iIndex = newIndex[(*conn).iIndex];				
				c.dValue = connValue;
				
				con.push_back(c);
				sumInterpolatory += connValue;
			}	
			
			
			double alpha = - (sumNeighbors / sumInterpolatory) / mnorm(A.getDiag(i));
			for(int j=0; j<con.size(); j++)
				con[j].dValue *= alpha;
			
			ASSERT2(con.size() > 0, "0 connections in point i = " << i << " ?");
			// connections hinzufügen
			P.setMatrixRow(i, &con[0], con.size());	
		}	
		else
		{
			ASSERT2(aggressiveCoarsening != 0, "no aggressive Coarsening but node " << i << " is fine and indirect??");
		}		
	}	
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateIndirectProlongation:
//-------------------------
//! Assume Prolongation of all normal fine nodes is already computed, it calculates the Interpolation of
//! fineIndirect nodes with Matrix_type A and Coarse/FineIndirect markers in grid[i].isCoarse/isFineIndirect
//!
//! Probably this is not the fastest way to do this:
//! One could create the graph 1 directly with indirect interpolation, then coarse, and then
//! calc interpolation. For fine nodes with no coarse neighbors, calc fine neighbors' interpolation,
//! then calc indirect interpolation. check if interpolation already calculated by looking at P.iNrOfConnections[i].
//!
//! uses nodeinfo *grid.
//!	@param P				Matrix P: here goes the calculated prolongation
//! @param	A				Matrix A: matrix for which to calculate prolongation on next level
//! @param	newIndex		newIndex of coarse Node i on next coarser level
//! @param	iNrOfCoarse		nr of coarse nodes on this level
//! @param posInConnections		array of size A.getLength() for speedup of neighbor-neighbor-calculation inited with -1.
template<typename Matrix_type, typename Vector_type>
void amg<Matrix_type, Vector_type>::CreateIndirectProlongation(SparseMatrix<double> &P, const Matrix_type &A, int *newIndex, int *posInConnections)
{
	ASSERT2(aggressiveCoarsening, "indirect interpolation only for aggressive coarsening");
	vector<SparseMatrix<double>::connection > con, con2;
	vector<int> nrOfPaths;
	con.reserve(255); con2.reserve(255); nrOfPaths.reserve(255);
	SparseMatrix<double>::connection c;
	//P.print();
	// INDIRECT INTERPOLATION
	
	for(int i=0; i<A.getLength(); i++)
	{
		if(!grid[i].isFineIndirect() || A[i].isUnconnected())
			continue;
		// a fine node with INDIRECT INTERPOLATION
		
		// calculate min offdiag-entry
		double dmax = 0;
		typename Matrix_type::cRowIterator conn = A.beginRow(i); ++conn;
		for(; !conn.isEnd(); ++conn)
			if(dmax < mnorm((*conn).dValue)) dmax = mnorm((*conn).dValue);	
		
		con.clear();
		con2.clear();
		nrOfPaths.clear();
		
		double sumInterpolatory=0, sumNeighbors=0;
		
		conn.rewind(); ++conn; // skip diagonal
		for(; !conn.isEnd(); ++conn)
		{
			double connValue = -mnorm((*conn).dValue);
			//sumNeighbors += connValue;
			if(-connValue < theta * dmax)
				continue;			
			
			int indexN = (*conn).iIndex;
			
			if(grid[indexN].isCoarse())
			{
				// add to interpolatory set
				c.iIndex = newIndex[(*conn).iIndex];
				c.dValue = connValue;
				con.push_back(c);
				sumInterpolatory += connValue;
				sumNeighbors += connValue;
				// continue oder abbruch???
				continue;
			}
			
			ASSERT2(grid[indexN].isFineIndirect() == FALSE, "indirect fine index " << i << " cannot have indirect fine neighbors " << indexN << "!");
			
			typename SparseMatrix<double>::rowIterator conn2 = P.beginRow(indexN); // !!! P
			++conn2; // skip diag TODO: true?
			for(; !conn2.isEnd(); ++conn2)
			{					
				int indexNN = (*conn2).iIndex;							
				int pos = posInConnections[indexNN];
				if(pos == -1)
				{
					posInConnections[indexNN] = con2.size();
					c.iIndex = indexNN;
					assign_mult(c.dValue, connValue, (*conn2).dValue);
					con2.push_back(c);
					nrOfPaths.push_back(1);
				}
				else
				{
					add_mult(con2[pos].dValue, connValue, (*conn2).dValue);
					nrOfPaths[pos]++;
				}						
			}
		}
		
		
		for(int j=0; j<con2.size(); j++)
		{			
			if(nrOfPaths[j] >= aggressiveCoarseningNrOfPaths)
			{
				con.push_back(con2[j]);
				sumInterpolatory += con2[j].dValue;
			}
			sumNeighbors += con2[j].dValue; // ???
			
			// reset posInConnections
			posInConnections[con2[j].iIndex] = -1;
		}
		
		ASSERT2(sumInterpolatory != 0.0, " numerical unstable?");
		
		double alpha = - (sumNeighbors / sumInterpolatory) / mnorm(A.getDiag(i));
		for(int j=0; j<con.size(); j++)
			con[j].dValue *= alpha;
		
		// connections hinzufügen
		P.setMatrixRow(i, &con[0], con.size());
		
	}
	P.finish();	
	//P.print();
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
void amg<Matrix_type, Vector_type>::createAMGLevel(Matrix_type &AH, SparseMatrix<double> &R, const Matrix_type &A, SparseMatrix<double> &P, int level)
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
	
	//cout << "AMG starts." << endl;
	// construct strong couplings-graph
	//cout << "construct strong couplings-graph" << endl;
	
	maxheap<nodeinfo> PQ(A.getLength(), grid);
	
	std::vector<bool>	coarse(A.getLength());
	int unassigned = A.getLength();
	
	int *newIndex = new int[A.getLength()];
	memset(newIndex, -1, sizeof(int)*A.getLength());
	
	int *posInConnections = new int[A.getLength()];
	memset(posInConnections, -1, sizeof(int)*A.getLength());
	cgraph graph(A.getLength(), A.getTotalNrOfConnections());
	
	/////////////////////////////////////////
	
	cout << "building graph... "; cout.flush();
	if(bTiming) SW.start();
	// CreateGraph(C, graph, PQ, unassigned);
	CreateGraph(A, graph, PQ, unassigned);	
	if(bTiming) SW.printTimeDiff();
	
	/* for(int i=0; i<A.getLength(); i++)
	 {
	 cout << i << " (" << GetPosForIndexAtLevel(i, level).x << " " << GetPosForIndexAtLevel(i, level).y << ") ";
	 grid[i].print();
	 }//*/
	//PQ.print();
	
	
	cout << "coarsening... "; cout.flush();	
	
	if(bTiming) SW.start();
	int iNrOfCoarse = Coarsen(graph, PQ, newIndex, unassigned, false, A);	
	
	if(bTiming) { SW.printTimeDiff();}
	
	cout << iNrOfCoarse << " Coarse Nodes. ";
	if(bTiming) cout << endl;
	
	
	/////////////////////////////////////////
	
	if(aggressiveCoarsening && level == 0)
	{
		cgraph graph2(A.getLength(), A.getTotalNrOfConnections());
		
		cout << "building graph2... "; cout.flush();
		if(bTiming) SW.start();
		
		CreateGraph2(graph, graph2, PQ, unassigned, posInConnections);	
		iNrOfCoarse -= unassigned;
		if(bTiming) SW.printTimeDiff();
		
		/*for(int i=0; i<A.getLength(); i++)
		 {
		 cout << i << " (" << GetPosForIndexAtLevel(i, level).x << " " << GetPosForIndexAtLevel(i, level).y << ") " ;
		 grid[i].print();
		 }//*/
		//PQ.print();
		
		
		cout << "coarsening2... "; cout.flush();	
		
		if(bTiming) SW.start();
		int newCoarseNodes = Coarsen(graph2, PQ, newIndex, unassigned, true, A);	
		iNrOfCoarse += newCoarseNodes;
		
		if(bTiming) { SW.printTimeDiff();}
		
		cout << iNrOfCoarse << " Coarse Nodes.";
		if(bTiming) cout << endl;
	}
	
	/////////////////////////////////////////
	
	vec1[level+1] = new Vector_type (iNrOfCoarse, "AMG:tempvec 1");
	vec1[level+1]->level = level+1;
	vec2[level+1] = new Vector_type (iNrOfCoarse, "AMG:tempvec 2");
	vec2[level+1]->level = level+1;
	cout << "created vec1 on level" << level +1 << endl;
	
	
	for(int i=0, c=0; i<A.getLength(); i++)
		if(grid[i].isCoarse())
		{			
			setSize((*vec1[level+1])[c], getRows(A[i][0].dValue));
			setSize((*vec2[level+1])[c], getRows(A[i][0].dValue));
			c++;
		}
	
	
	parentIndex[level+1] = new int[iNrOfCoarse];
	for(int i=0; i<A.getLength(); i++)
		if(grid[i].isCoarse())
			parentIndex[level+1][ newIndex[i] ] = i;
	
	/*for(int i=0; i<A.getLength(); i++)
	 {
	 cout << i << " (" << GetPosForIndexAtLevel(i, level).x << " " << GetPosForIndexAtLevel(i, level).y << ") ";
	 grid[i].print();
	 }//*/
	
#ifdef AMG_PRINT_COARSENING
	printCoarsening(level);
#endif
	
	
	// construct prolongation P = I_{2h->h}
	// ASSUME THAT Matrix_type IS SYMMETRIC!!!
	
	/////////////////////////////////////////
	
	cout << "prolongation... "; cout.flush();
	
	P.fromlevel = level+1;
	P.tolevel = level;
	
	if(bTiming) SW.start();
	CreateProlongation(P, A, newIndex, iNrOfCoarse);
	if(aggressiveCoarsening)
		CreateIndirectProlongation(P, A, newIndex, posInConnections);
	if(bTiming) SW.printTimeDiff();
	
#ifdef AMG_PRINT_P
	P.print();
#endif
	
	cout << "restriction... "; cout.flush();
	if(bTiming) SW.start();
	// construct restriction R = I_{h -> 2h}		
	R.createAsTransposeOf(P); // already finished
	R.name = "AMG:R";
	//R.print("R");
	if(bTiming) SW.printTimeDiff();	
	
#ifdef AMG_PRINT_R
	R.print();
#endif
	
	/////////////////////////////////////////	
	// create Galerkin product
	
	cout << "galerkin product... "; cout.flush();
	if(bTiming) SW.start();
	
	// AH = R A P
	AH.createAsMultiplyOf(R, A, P, posInConnections);
	//createGalerkinMatrix(AH, R, A, P, posInConnections);
	
	/*	SparseMatrix<double> RA;
	 RA.createAsMultiplyOf(R, A);
	 AH.createAsMultiplyOf(RA, P); */
	
	/*	SparseMatrix<double> AP;
	 AP.createAsMultiplyOf(A, P);
	 AH.createAsMultiplyOf(R, AP);*/
	
	
	if(bTiming) SW.printTimeDiff();
	AH.name = "AMG:A";
	AH.fromlevel = level+1;
	AH.tolevel = level+1;
	AH.finish();
	
#ifdef AMG_PRINT_AH
	AH.print();
#endif
	
	
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
	
	
#ifdef AMG_WRITE_MATRICES_PATH
	if(this->A[0]->getLength() < AMG_WRITE_MATRICES_MAX)
	{
		cout << "write matrices";
		P.writeToFile((string(AMG_WRITE_MATRICES_PATH) + "P" + nrstring(level) + ".mat").c_str()); cout << "."; cout.flush();
		R.writeToFile((string(AMG_WRITE_MATRICES_PATH) + "R" + nrstring(level) + ".mat").c_str()); cout << "."; cout.flush();
		AH.writeToFile((string(AMG_WRITE_MATRICES_PATH) + "A" + nrstring(level+1) + ".mat").c_str()); cout << "."; cout.flush();
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
// amg<Matrix_type, Vector_type>::onlyOneLevel
//----------------
//! for testing. creates only one AMG level without direct solvers
//! @param A	matrix A.
template<typename Matrix_type, typename Vector_type>
bool amg<Matrix_type, Vector_type>::onlyOneLevel(const Matrix_type& A_)
{
	used_levels = 2;
	const Matrix_type *pA = &A_;
	grid = new nodeinfo[pA->getLength()];
	A[0] = const_cast<Matrix_type*> (pA);
	
	int i=0;
	A[i+1] = new Matrix_type();
	createAMGLevel(*A[i+1], R[i], *A[i], P[i], i);
	
	//	A[1]->print();
	return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// amg<Matrix_type, Vector_type>::init
//----------------
//! creates MG Hierachy for with Matrix_type A and temporary vectors for higher levels
//! @param A	matrix A.
template<typename Matrix_type, typename Vector_type>
bool amg<Matrix_type, Vector_type>::init(const Matrix_type& A_)
{
	cout << "Starting AMG Setup." << endl << endl;
	stopwatch SW;
	const Matrix_type *pA = &A_;
	grid = new nodeinfo[pA->getLength()];
	A[0] = const_cast<Matrix_type*> (pA);
	
	A[0]->fromlevel = 0;
	A[0]->tolevel = 0;
	
	
#ifdef AMG_WRITE_MATRICES_PATH	
	if(A[0]->getLength() < AMG_WRITE_MATRICES_MAX)
	{
		cout << "write matrix A...";
		A[0]->writeToFile((string(AMG_WRITE_MATRICES_PATH) + "A" + nrstring(0) + ".mat").c_str());
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
	double totallength;
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
amg<Matrix_type, Vector_type>::amg()
{
	used_levels = 0;
	max_levels = 10;
	aggressiveCoarsening = 0;
	aggressiveCoarseningNrOfPaths = 2; // A2
	
	nu1 = 2;
	nu2 = 2;
	gamma = 1;
	
	//	if(never_happens) printCoarsening(0,0);
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
	}
}

#ifndef DEBUG

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MGCycle:
//-------------------------
//! calculates one AMG cycle
//! @param		x		
//! @param	b	
//! @param	level	
template<typename Matrix_type, typename Vector_type>
double amg<Matrix_type, Vector_type>::MGCycle(Vector_type &x, const Vector_type &b, int level)
{	
	ASSERT2(x.getLength() == b.getLength() && b.getLength() == A[level]->getLength(), 
			"x (" << x << "), b (" << b << "), or A (" << *A[level] << ") have different length");
	
	/*spaceout(level);
	 cout << "AMG MGCycle, level " << level;
	 cout.flush(); //*/
	const Matrix_type &Ah = *(A[level]);
	
	if(level == used_levels-1)
	{
		coarseSolver.solve(b, x);
		return 0.1e-14;
	}
	
	//cout << level << " pre presmooth " << norm(b - Ah*x) << endl;
	
	for(int i=0; i < nu1; i++)
		smoother[level].iterate(x, b);
	Vector_type &r = *vec3[level];
	
	r = b - Ah*x;
	//cout << level << " presmooth " << norm(r) << endl;
	
	Vector_type &rH = *vec1[level+1];
	Vector_type &eH = *vec2[level+1];
	
	rH = R[level]*r; //restriction(r, level);
	//P[level].transposeMult(rH, r);
	
	eH = 0.0;
	
	if(level+1 == used_levels-1)
		MGCycle(eH, rH, level+1);
	else
		for(int i=0; i<gamma; i++)
			MGCycle(eH, rH, level+1);
	
	/*	Matrix_type &aa = (*A[level+1]);
	 Vector_type y(x.getLength()), rr(x.getLength());
	 y = P[level]*eH;
	 x += P[level]*eH; //interpolate(eH, level);	
	 rr = b-Ah*x;	
	 cout << level << " pre postsmooth " << norm(b - Ah*x) << endl;*/
	
	x += P[level]*eH; //interpolate(eH, level);	
	
	double res;
	if(level != 0)
		for(int i=0; i < nu2; i++)
			res = smoother[level].iterate(x, b);
	
	//cout << level << " postsmooth " << norm(b - Ah*x) << endl;
	
	return res;
}
#else

template<typename Matrix_type, typename Vector_type>
double amg<Matrix_type, Vector_type>::MGCycle(Vector_type &x, const Vector_type &b, int level)
{	
	ASSERT2(x.getLength() == b.getLength() && b.getLength() == A[level]->getLength(), 
			"x (" << x << "), b (" << b << "), or A (" << *A[level] << ") have different length");
	
	
	//spaceout(2*level); cout << "AMG MGCycle, level " << level << endl;
	cout.flush(); 
	const Matrix_type &Ah = *(A[level]);
	
	if(level == used_levels-1)
	{
		coarseSolver.solve(b, x);
		return 0.1e-14;
	}
	
	double pre1 = norm(b - Ah*x);
	
	for(int i=0; i < nu1; i++)
		smoother[level].iterate(x, b);
	Vector_type &r = *vec3[level];
	
	r = b - Ah*x;
	double presmoothreduction = norm(r)/pre1;
	//spaceout(2*level); cout << "presmooth reduction: " << norm(r)/pre1 << endl;
	
	Vector_type &rH = *vec1[level+1];
	Vector_type &eH = *vec2[level+1];
	
	rH = R[level]*r; 
	
	eH = 0.0;
	
	if(level+1 == used_levels-1)
		MGCycle(eH, rH, level+1);
	else
		for(int i=0; i<gamma; i++)
			MGCycle(eH, rH, level+1);
	
	
	x += P[level]*eH; //interpolate(eH, level);	
	double pre2 = norm(b - Ah*x);
	
	double res;
	//if(level != 0)
		for(int i=0; i < nu2; i++)
			res = smoother[level].iterate(x, b);
	
	double post = norm(b - Ah*x);
	
	
	//spaceout(2*level); cout << "postsmooth reduction: " << post/pre2 << endl;
	//spaceout(2*level); cout << "level reduction: " << post/pre1 << endl;
	 cout << "	[" << level << "] prered: " << presmoothreduction << " postred: " << post/pre2
		<< " levelred: " << post/pre1 << endl;
	
	return res;
}


#endif