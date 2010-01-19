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


#if 0
#define AMG_PRINT_COARSEN
#define AMG_PRINT_COARSENING

#define AMG_PRINT_P
#define AMG_PRINT_R
#define AMG_PRINT_AH
#define AMG_PRINT_GRAPH
#endif
//*/
// writeMatrices:
//----------------
// writes A to pathAndName + "A" + level + ".mat for all levels, also R
template<typename entry_type, typename Vector_type>
void amg<entry_type, Vector_type>::writeMatrices(const char *pathAndName)
{
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

template<typename entry_type, typename Vector_type>
void amg<entry_type, Vector_type>::restriction(Vector_type *pto, const Vector_type &from, int fromlevel) // h -> 2h
{	
	Vector_type &to = *pto;
	ASSERT2(to.getLength() == R[fromlevel].getLength() && from.getLength() == P[fromlevel].getLength(), "lenghts not matching");
	
	to = R[fromlevel] * from;
}

template<typename entry_type, typename Vector_type>
void amg<entry_type, Vector_type>::interpolate(Vector_type *pto, const Vector_type &from, int tolevel) // 2h -> h prolongate
{
	Vector_type &to = *pto;
	ASSERT2( to.getLength() == P[tolevel].getLength() && from.getLength() == R[tolevel].getLength(), "lengths not matching");
	
	to = P[tolevel]*from;
}

// printCoarsening:
//----------------
// Debug output. Writes position of Coarse nodes in coarse<level>.dat, and fine in fine<level>.dat for display in gnuplot
template<typename entry_type, typename Vector_type>
void amg<entry_type, Vector_type>::printCoarsening(int level, int n)
{
  
	fstream fcoarse((string("/Users/mrupp/matrices/coarse") + nrstring(level) + ".dat").c_str(), ios::out);	
	fstream ffine  ((string("/Users/mrupp/matrices/fine") + nrstring(level) + ".dat").c_str(), ios::out);	
	//fstream funass ((string("unass") + nrstring(level) + ".dat").c_str(), ios::out);	
	
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
// CreateGraph:
//-------------------------
// Create graph of strong connections from A, calculate ratings in grid[i].rating,
// build up priority queue PQ. unassigned: nr of nodes to be assigned by coarsening algorihm
template<typename entry_type, typename Vector_type> // template<typename conn_matrix> // const conn_matrix &C
void amg<entry_type, Vector_type>::CreateGraph(const matrix_type &A, cgraph &graph, maxheap<nodeinfo> &PQ, int &unassigned)
{
	unassigned = 0;
	for(int i=0; i< A.getLength(); i++)
	{
		graph.init(i);
		if(A[i].isUnconnected())
			continue;
		
		double dmax = 0;
		typename matrix_type::cRowIterator conn = A.beginRow(i); ++conn;
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
	// and NOT all nodes which interpolate the coarse node
	graph.transpose();
	for(int i=0; i < A.getLength(); i++)
	{
		if(A[i].isUnconnected())
			grid[i].setFine();
		else
		{
			//ASSERT2(graph.iNrOfConnections[i] > 0, "node " << i << " has " << graph.iNrOfConnections[i] << " connections?");
			grid[i].rating = graph.iNrOfConnections[i];
			PQ.insertItem(i);
			unassigned++;
		}
	}	
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateGraph2:
//-------------------------
// Create the graph2, which consists of coarse nodes in graph
// two coarse nodes a and b in graph2 are connected, if
// - there exist at least 2 ways of length 2 from a to b (A2-Coarsening)
// - there exist a ways of length 2 from a to b (A1-Coarsening)
// because of the coarsening process, every coarse node has only fine neighbors in graph,
// that means those ways are from a coarse node over a fine node to a coarse node.
template<typename entry_type, typename Vector_type>
void amg<entry_type, Vector_type>::CreateGraph2(cgraph &graph, cgraph &graph2, maxheap<nodeinfo> &PQ, int &unassigned, int *posInConnections)
{	
	vector<int> connection(255);
	vector<int> nrOfPaths(255);
	
	PQ.reset();
	//graph.print();
	unassigned=0;
	for(int b=0; b < graph.size; b++)
	{
		graph2.init(b);
		if(grid[b].isFine())
			continue;
		
		connection.clear(); nrOfPaths.clear();
		// first calculate all nodes reachable with paths of length 2
		
		// ! i is coarse -> has only fine neighbors
		for(int i=0; i<graph.iNrOfConnections[b]; i++)
		{
			int indexN = graph.conn[b][i];					
			for(int j=0; j<graph.iNrOfConnections[indexN]; j++)
			{
				int indexNN = graph.conn[indexN][j];
				
				if(grid[indexNN].isFine() || indexNN == b)
					continue;
				int pos = posInConnections[indexNN];
				if(pos == -1)
				{					
					pos = posInConnections[indexNN]= connection.size();					
					connection.push_back(indexNN);
					nrOfPaths.push_back(1);					
				}
				else
					nrOfPaths[pos]++;
			}
		}
		
		// then sort out those which were reached at #aggressiveCoarseningNrOfPaths (2 or 1)
		grid[b].rating = 0;
		for(int i=0; i<connection.size(); i++)
		{
			if(nrOfPaths[i] >= aggressiveCoarseningNrOfPaths)
			{
				graph2.setConnection(b, connection[i]);
				grid[b].rating ++;
			}
			// reset posInConnections for further use
			posInConnections[connection[i]] = -1;
		}
		
		// add node with rating to priority queue
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
// Coarsens the graph with ratings of nodes in grid[i].rating, set up in a priority queue PQ
// newIndex[i]: store here new index of this node in coarse grid (>0, if fine < 0)
// int unassigned: nr of nodes to assign
// bool bIndirect: if true, this is 2nd stage of Aggressive Coarsening, then fine nodes get marker "IndirectFine"
//				   instead of just "fine". Used later in CreateProlongation and CreateIndirectProlongation
template<typename entry_type, typename Vector_type>
int amg<entry_type, Vector_type>::Coarsen(cgraph &graph, maxheap<nodeinfo> &PQ, int *newIndex, int unassigned, bool bIndirect, const matrix_type &A)
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
		int nrOfConnections = graph.iNrOfConnections[best];
		
		for(int i=0; i < nrOfConnections; i++)
		{
			//cout << graph.conn[best][i] << " ";
			int indexN = graph.conn[best][i];
			if(grid[indexN].isAssigned())
				continue;					
			PQ.remove(indexN);			
		}
		
		// mark neighbors as fine
		//cout << " fine: ";
		
		
		for(int i=0; i<nrOfConnections; i++)
		{
			int indexN = graph.conn[best][i];
			
			
#ifdef AMG_PRINT_COARSEN
			cout << indexN << " [" << GetOriginalIndex(A.tolevel, indexN) << "] "<< " ";
			if(grid[indexN].isAssigned())
				cout << grid[indexN].isCoarse()) ? "(c) " : "(f) ";
#endif
			
			if(grid[indexN].isAssigned())
				continue;
			
			if(bIndirect) grid[indexN].setFineIndirect();
			else grid[indexN].setFine();
			
			

			unassigned--;
			
			// increase rating of neighbors of neighbors
			int nrOfConnectionsN = graph.iNrOfConnections[indexN];
			for(int j=nrOfConnectionsN-1; j>0; j--)
			{
				int indexNN = graph.conn[indexN][j];
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
	return iNrOfCoarse;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateProlongation:
//-------------------------
// Calculates Prolongatin P with matrix_type A and coarse/fine markers in grid[i].isFine/isCoarse
template<typename entry_type, typename Vector_type>
void amg<entry_type, Vector_type>::CreateProlongation(SparseMatrix<double> &P, const matrix_type &A, int *newIndex, int iNrOfCoarse)
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
			
			// calc min off-diag-entry
			double dmax = 0;
			typename matrix_type::cRowIterator conn = A.beginRow(i); ++conn;
			for(; !conn.isEnd(); ++conn)
				if(dmax < mnorm((*conn).dValue)) dmax = mnorm((*conn).dValue);			
			
			// calculate interpolation coefficients
			double sumNeighbors =0, sumInterpolatory=0;
			
			con.clear(); double connValue;
			
			conn.rewind(); ++conn; // skip diagonal
			for(; !conn.isEnd(); ++conn)
			{
				connValue = -mnorm((*conn).dValue);
				sumNeighbors += connValue;	
				
				if(-connValue < theta * dmax || grid[(*conn).iIndex].isFine())
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
// Assume Prolongation of all normal fine nodes is already computed, it calculates the Interpolation of
// fineIndirect nodes with matrix_type A and Coarse/FineIndirect markers in grid[i].isCoarse/isFineIndirect
//
// Probably this is not the fastest way to do this:
// One could create the graph 1 directly with indirect interpolation, then coarse, and then
// calc interpolation. For fine nodes with no coarse neighbors, calc fine neighbors' interpolation,
// then calc indirect interpolation. check if interpolation already calculated by looking at P.iNrOfConnections[i].
template<typename entry_type, typename Vector_type>
void amg<entry_type, Vector_type>::CreateIndirectProlongation(SparseMatrix<double> &P, const matrix_type &A, int *newIndex, int *posInConnections)
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
		typename matrix_type::cRowIterator conn = A.beginRow(i); ++conn;
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
			++conn2; // skip diag
			for(; !conn2.isEnd(); ++conn2)
			{					
				int indexNN = (*conn2).iIndex;							
				int pos = posInConnections[indexNN];
				if(pos == -1)
				{
					posInConnections[indexNN] = con2.size();
					c.iIndex = indexNN;
					c.dValue = connValue * (*conn2).dValue;
					con2.push_back(c);
					nrOfPaths.push_back(1);
				}
				else
				{
					con2[pos].dValue += connValue * (*conn2).dValue;
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
// createGalerkinMatrix:
//-------------------------
// Calculates AH = R A P. posInConnections only needed for speedup (has to be -1 forall i).
template<typename entry_type, typename Vector_type>
void amg<entry_type, Vector_type>::createGalerkinMatrix(matrix_type &AH, const SparseMatrix<double> &R, const matrix_type &A, const SparseMatrix<double> &P, int *posInConnections)
{
	//cout << endl << " Creating Galerkin matrix_type..." << endl;
	int n = R.getLength();
	
	// speedup with array posInConnections, needs n memory (also used in CreateGraph2 and CreateIndirectProlongation).
	// 1000x1000 ninepoint: old version: 2444 ms, new 1200 ms, 1000 ms, 900 ms, 950 ms
	
	// posInConnections[i]: index in the connections for current row.
	// has to be -1 for all nodes	
	
	vector<typename matrix_type::connection > con(255);
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
		typename matrix_type::connection c;
		c.iIndex = i;
		c.dValue = 0.0;
		con.push_back(c);
		
		for(typename SparseMatrix<double>::cRowIterator itR(R, i); !itR.isEnd(); ++itR)
		{
			r = (*itR).dValue;
			if(r == 0.0) continue;
			
			for(typename matrix_type::cRowIterator itA(A, (*itR).iIndex); !itA.isEnd(); ++itA)
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
// createAMGLevel:
//-------------------------
// 
template<typename entry_type, typename Vector_type>
void amg<entry_type, Vector_type>::createAMGLevel(matrix_type &AH, SparseMatrix<double> &R, const matrix_type &A, SparseMatrix<double> &P, int level)
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
	printCoarsening(level, A.getLength());
#endif
	
	
	// construct prolongation P = I_{2h->h}
	// ASSUME THAT matrix_type IS SYMMETRIC!!!
	
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
	createGalerkinMatrix(AH, R, A, P, posInConnections);
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
	cout << "	nnz: " << nnz << " n2: " << (double(AH.getLength())*double(AH.getLength()))
	<< " Density: " << double(nnz)/(double(AH.getLength())*double(AH.getLength()))*100.0 << "% nnz/n: " << nnz/(double)AH.getLength() << endl;
	
	cout << " level "; SWwhole.printTimeDiff();  cout << endl;
	cout << endl;
	cout.flush();
	
	cout << A << endl;
	cout << P << endl;
	cout << R << endl;
	
	//P.print();
	
	cout << AH << endl;
	cout << endl;//*/
	
	delete[] posInConnections;
	delete[] newIndex;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// amg<entry_type, Vector_type>::onlyOneLevel
//----------------
// for testing
template<typename entry_type, typename Vector_type>
bool amg<entry_type, Vector_type>::onlyOneLevel(const matrix_type& A_)
{
	used_levels = 2;
	const matrix_type *pA = &A_;
	grid = new nodeinfo[pA->getLength()];
	A[0] = const_cast<matrix_type*> (pA);
	
	int i=0;
	A[i+1] = new matrix_type();
	createAMGLevel(*A[i+1], R[i], *A[i], P[i], i);

//	A[1]->print();
	return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// amg<entry_type, Vector_type>::init
//----------------
// creates MG Hierachy for with matrix_type A and temporary vectors for higher levels
template<typename entry_type, typename Vector_type>
bool amg<entry_type, Vector_type>::init(const matrix_type& A_)
{
	stopwatch SW;
	const matrix_type *pA = &A_;
	grid = new nodeinfo[pA->getLength()];
	A[0] = const_cast<matrix_type*> (pA);
	
	A[0]->fromlevel = 0;
	A[0]->tolevel = 0;
	
	
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
		
		A[i+1] = new matrix_type();
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
	used_levels = i+1;
	SW2.printTimeDiff();
	cout << endl << endl;
	
	cout << "AMG Setup finished. Used Levels: " << used_levels << ". ";
	SW.printTimeDiff(); cout << endl << endl;
	delete [] grid;
	
	return true;
}

template<typename entry_type, typename Vector_type>
amg<entry_type, Vector_type>::amg()
{
	used_levels = 0;
	max_levels = 10;
	aggressiveCoarsening = 0;
	aggressiveCoarseningNrOfPaths = 2; // A2
	
	nu1 = 2;
	nu2 = 2;
	gamma = 1;
	
	if(never_happens) printCoarsening(0,0);
}

template<typename entry_type, typename Vector_type>
amg<entry_type, Vector_type>::~amg()
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
// 
template<typename entry_type, typename Vector_type>
double amg<entry_type, Vector_type>::MGCycle(Vector_type &x, const Vector_type &b, int level)
{
	
	ASSERT2(x.getLength() == b.getLength() && b.getLength() == A[level]->getLength(), 
			"x (" << x << "), b (" << b << "), or A (" << *A[level] << ") have different length");
	
	/*spaceout(level);
	 cout << "AMG MGCycle, level " << level;
	 cout.flush(); //*/
	const matrix_type &Ah = *(A[level]);
	
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
	
/*	matrix_type &aa = (*A[level+1]);
	Vector_type y(x.getLength()), rr(x.getLength());
	y = P[level]*eH;
	x += P[level]*eH; //interpolate(eH, level);	
	rr = b-Ah*x;	
	cout << level << " pre postsmooth " << norm(b - Ah*x) << endl;*/
	
	x += P[level]*eH; //interpolate(eH, level);	
	
	double res;
	for(int i=0; i < nu2; i++)
		res = smoother[level].iterate(x, b);
	
	//cout << level << " postsmooth " << norm(b - Ah*x) << endl;
	
	return res;
}

