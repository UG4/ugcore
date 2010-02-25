/*
 *  amg.hpp
 *  flexamg
 *
 *  Created by Martin Rupp on 03.12.09.
 *  Copyright 2009 . All rights reserved.
 *
 */
#pragma once

extern int *parentIndex[AMG_MAX_LEVELS];
//#define GRAPH_WITH_LOCAL_INVERSE


#define AMG_WRITE_MATRICES_PATH "/Users/mrupp/matrices/AMG_"
#define AMG_WRITE_MATRICES_MAX (200*200)


//#define AMG_PRINT_INDIRECT

//#define AMG_PRINT_GRAPH

#define AMG_WRITE_COARSENING
#define AMG_WRITE_GRAPH

#if 0

#define AMG_PRINT_COARSENING
#define AMG_PRINT_P
#define AMG_PRINT_R
#define AMG_PRINT_AH

#define AMG_PRINT_COARSEN_RATINGS
#define AMG_PRINT_COARSEN

#endif

inline double amg_value(double d)
{
	return d;
}

template<typename T>
double amg_value(T &d)
{
	return -d.norm();
}
//
// writeMatrices:
//----------------
//! writes A to pathAndName + "A" + level + ".mat for all levels, also R
//! @param pathAndName	path with name of matrix. for example "/Users/username/matrices/mat"
template<typename Matrix_type, typename Vector_type>
void amg<Matrix_type, Vector_type>::writeMatrices(const char *pathAndName)
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
#ifndef GRAPH_WITH_LOCAL_INVERSE
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
#else
	
	fstream cout;
	
	vector<int> indices;
	blockDenseMatrix<double, variableStorage> SM(indices.size(), indices.size());

	unassigned = 0;
	for(int i=0; i< A.getLength(); i++)
	{
		graph.init(i);
		if(A[i].isUnconnected())
			continue;
		
		A.getNeighborhood(i, 2, indices);
		SM.setSize(indices.size(), indices.size());
		cout << '\t' << i << endl;
		A.get(SM, indices, indices);
		
		
	/*	for(int c=0; c<indices.size(); c++)
			cout << indices[c] << " ";
		cout << endl;*/
		
		SM.invert();
		
		double dmax = 0;
		int ind = find(indices.begin(), indices.end(), i) - indices.begin();
		
		cout << endl;
		double alpha = 1/mnorm(SM(ind, ind));
		for(int c=0; c<SM.getCols(); c++)
		{
			SM(ind,c)*=alpha;
			if(c!=ind && dmax < SM(ind,c))
				dmax = SM(ind,c); 
			cout << setw(13) << indices[c];
		}
		cout << endl;
		//SM.print();
		
		for(int c=0; c<SM.getCols(); c++)
		{
			
			if(c!= ind && 0.4*dmax < SM(ind,c))
			{
				graph.setConnection(i, indices[c]);	
				cout << c << " ";
				cout << "*";
			}
			else cout << " ";
			cout << setw(12) << left << SM(ind, c);
		}
		cout << endl << "dmax is " << dmax << endl << endl;
		aggressiveCoarsening = true;
	}

	
#endif
#ifdef AMG_PRINT_GRAPH
	graph.print();
#endif

#ifdef AMG_WRITE_GRAPH
	graph.writeToFile((string(AMG_WRITE_MATRICES_PATH) + "G" + nrstring(A.tolevel) + ".mat").c_str(), A.tolevel);
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
void amg<Matrix_type, Vector_type>::CreateGraph2(cgraph &graph, cgraph &graph2, maxheap<nodeinfo> &PQ, int &unassigned, int &iNrOfCoarse, int *posInConnections, int *newIndex)
{		
	vector<int> connection(255);
	vector<int> nrOfPaths(255);
	
	nodeinfo *grid2 = new nodeinfo[graph.getLength()];
	
	iNrOfCoarse = 0;
	PQ.create(graph.getLength(), grid2);
	//graph.print();
	unassigned=0;
	for(int b=0; b < graph.getLength(); b++)
	{
		graph2.init(b);
		if(grid[b].isFineDirect())
		{
			grid2[b].setFineDirect();
			continue;
		}
		
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
				
				if(indexNN == b || grid[indexNN].isFineDirect())
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
		grid2[b].rating = 0;
		for(int i=0; i<connection.size(); i++)
		{
			if(nrOfPaths[i] >= aggressiveCoarseningNrOfPaths)
			{
				// add connection b -> node
				graph2.setConnection(b, connection[i]);
				// increase rating of b
				grid2[b].rating ++;
			}
			// reset posInConnections for further use
			posInConnections[connection[i]] = -1;
		}
		
		// add node with rating > 0 to priority queue
		if(grid2[b].rating > 0)
		{
			PQ.insertItem(b);
			unassigned++;
		}
		else
		{
			grid2[b].setCoarse();
			newIndex[b] = iNrOfCoarse++;
		}		
	}	

	delete[] grid;
	grid = grid2;

	
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
int amg<Matrix_type, Vector_type>::Coarsen(cgraph &graph, maxheap<nodeinfo> &PQ, int *newIndex, int unassigned, int &iNrOfCoarse, const Matrix_type &A)
{
#ifdef AMG_WRITE_COARSENING
	fstream fstr((string(AMG_WRITE_MATRICES_PATH) + "A" + nrstring(A.tolevel) + ".mat").c_str(), ios::out|ios::app);
	fstr << "c" << endl;
#endif
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
#ifdef AMG_WRITE_COARSENING
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
			
#ifdef AMG_PRINT_COARSEN
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
	// seems to work but doesnt help with convergence rates on complicated geometries???	
	cgraph TG;
	TG.createAsTransposeOf(graph);
	
	int nrOfFFCoarseNodes=0;
	vector<bool> marks(graph.getLength());
	
	// seems to work but doesnt help with convergence rates on complicated geometries???	
	
#ifdef AMG_WRITE_COARSENING
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
#ifdef AMG_WRITE_COARSENING
				fstr << GetOriginalIndex(A.tolevel, i) << endl << endl << GetOriginalIndex(A.tolevel, it()) << endl;
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
void amg<Matrix_type, Vector_type>::CreateProlongation(SparseMatrix<double> &P, const Matrix_type &A, int *newIndex, int iNrOfCoarse, int &unassigned)
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
			
			// calc min off-diag-entry
			double dmax = 0, connValue, maxConnValue = 0;
			double sumNeighbors =0, sumInterpolatory=0;
			typename Matrix_type::cRowIterator conn = A.beginRow(i); ++conn; // skip diag
			
			double diag = amg_value(A.getDiag(i));
			
			for(; !conn.isEnd(); ++conn)
			{
				connValue = amg_value((*conn).dValue);
				
				if(connValue > 0)
				{
					diag += connValue;
					continue;
				}
				
				sumNeighbors += connValue;
				
				if(dmax > connValue) 
					dmax = connValue;
				if(grid[(*conn).iIndex].isCoarse() && maxConnValue > connValue)
					maxConnValue = connValue;
				
			}
			

			double barrier;
			if(eps_truncation_of_interpolation > 0)  // Ruge/Stuebe A.7.2.4 truncation of interpolation
				barrier = min(theta*dmax, eps_truncation_of_interpolation*maxConnValue);
			else
				barrier = theta*dmax;
			
			con.clear();			
			
			conn.rewind(); ++conn; // skip diagonal
			for(; !conn.isEnd(); ++conn)
			{
				if(!grid[(*conn).iIndex].isCoarse()) continue;
				
				connValue = amg_value((*conn).dValue);	
				if(connValue > barrier)
					continue;
				c.iIndex = newIndex[(*conn).iIndex];   ASSERT1(c.iIndex >= 0);			
				c.dValue = connValue;
				
				con.push_back(c);
				sumInterpolatory += connValue;
			}	

			// connections hinzufügen
			if(con.size() > 0)
			{
				double alpha = - (sumNeighbors / sumInterpolatory) / diag;
				for(int j=0; j<con.size(); j++)
					con[j].dValue *= alpha;
				
				//ASSERT2(con.size() > 0, "0 connections in point i = " << i << " ?");
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
void amg<Matrix_type, Vector_type>::CreateIndirectProlongation(SparseMatrix<double> &P, const Matrix_type &A, int *newIndex, int *posInConnections, int unassigned)
{
	ASSERT2(aggressiveCoarsening, "indirect interpolation only for aggressive coarsening");
	vector<SparseMatrix<double>::connection > con, con2;
	vector<int> nrOfPaths;
	con.reserve(255); con2.reserve(255); nrOfPaths.reserve(255);
	SparseMatrix<double>::connection c;
	//P.print();
	// INDIRECT INTERPOLATION
	
	int oldUnassigned = -1;
	int pass=2;
	while(unassigned)
	{
#ifdef AMG_PRINT_INDIRECT
		cout << endl;
#endif
		cout << "Pass " << pass << ": ";
		for(int i=0; i<A.getLength() && unassigned > 0; i++)
		{
			if(!grid[i].isUnassignedFineIndirect() || A[i].isUnconnected())
				continue;

			double diag = amg_value(A.getDiag(i));
			// calculate min offdiag-entry
			double dmax = 0;
			typename Matrix_type::cRowIterator conn = A.beginRow(i); ++conn;
			for(; !conn.isEnd(); ++conn)
			{
				double connValue = amg_value((*conn).dValue);
				if(connValue > 0)
				{
					diag += connValue;
					continue;
				}
				if(dmax > connValue) 
					dmax = connValue;	
			}
			
			con.clear();
			con2.clear();
			nrOfPaths.clear();
			
			double sumInterpolatory=0, sumNeighbors=0;
			
			//cout << "indirect interpolating node " << i << endl;
			
			conn.rewind(); ++conn; // skip diagonal
			for(; !conn.isEnd(); ++conn)
			{
				int indexN = (*conn).iIndex;
				if(grid[indexN].isFineIndirectLevel(pass))
					continue;
				
				double connValue = amg_value((*conn).dValue);
				sumNeighbors += connValue;
				if(connValue > theta * dmax)
					continue;
				
				ASSERT2(!grid[indexN].isCoarse(), "Node " << i << " indirect, but neighbor " <<  indexN << " coarse?");
				
				//ASSERT2(grid[indexN].isFineIndirect() == FALSE, "indirect fine index " << i << " cannot have indirect fine neighbors " << indexN << "!");
				
				typename SparseMatrix<double>::rowIterator conn2 = P.beginRow(indexN); // !!! P
				++conn2;
				for(; !conn2.isEnd(); ++conn2)
				{					
					int indexNN = (*conn2).iIndex;							
					int pos = posInConnections[indexNN];
					
					if(pos == -1)
					{
						pos = posInConnections[indexNN] = con2.size();
						c.iIndex = indexNN; ASSERT1(c.iIndex >= 0);
						
						assign_mult(c.dValue, connValue, (*conn2).dValue);
						con2.push_back(c);
						//nrOfPaths.push_back(1);
					}
					else
					{
						add_mult(con2[pos].dValue, connValue, (*conn2).dValue);
						//nrOfPaths[pos]++;
					}
				}
			}
					
			for(int j=0; j<con2.size(); j++)
			{			
				//if(nrOfPaths[j] >= aggressiveCoarseningNrOfPaths)
				{
					con.push_back(con2[j]);

					sumInterpolatory += con2[j].dValue;
				}
				//sumNeighbors += con2[j].dValue; // ???
				
				// reset posInConnections
				posInConnections[con2[j].iIndex] = -1;
			}
			
			if(con.size() == 0)
				continue;
			
			unassigned --;
			
			grid[i].setFineIndirectLevel(pass);
#ifdef AMG_PRINT_INDIRECT
			cout << i << " ";
#endif
			//cout << endl;
			
			ASSERT2(sumInterpolatory != 0.0, " numerical unstable?");			
			double alpha =  /*1/sumInterpolatory; */ - (sumNeighbors / sumInterpolatory)/diag;
			for(int j=0; j<con.size(); j++)
			{
				//cout << con[j].dValue << " - N:" << sumNeighbors << " I: " << sumInterpolatory << " alpha: " << alpha << ". " << con[j].dValue*alpha << " : " << A.getDiag(i) << endl;
				con[j].dValue *= alpha;
			}
			
			// connections hinzufügen
			P.setMatrixRow(i, &con[0], con.size());
			
		}
		
		ASSERT2(unassigned != oldUnassigned, "Pass " << pass << ": Indirect Interpolation hangs at " << unassigned << " unassigned nodes.")
			
#ifdef AMG_PRINT_INDIRECT
		cout << "calculated, ";
#endif		
		cout << unassigned << " left. ";
		pass++;
		oldUnassigned = unassigned;
		//break;
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
	
	maxheap<nodeinfo> PQ(A.getLength(), grid);
	
	std::vector<bool>	coarse(A.getLength());
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
	
#ifdef AMG_PRINT_COARSEN_RATINGS
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
	
	// agressive Coarsening
	/////////////////////////////////////////	
#ifndef	GRAPH_WITH_LOCAL_INVERSE
	if(aggressiveCoarsening && level == 0)
	{
		// build graph 2
		//------------------

		cgraph graph2(A.getLength(), A.getTotalNrOfConnections());
		
		cout << "building graph2... "; cout.flush();
		if(bTiming) SW.start();
		
		memset(newIndex, -1, sizeof(int)*A.getLength());
		
		CreateGraph2(graph, graph2, PQ, unassigned, iNrOfCoarse, posInConnections, newIndex);	
		if(bTiming) SW.printTimeDiff();
		
		/*for(int i=0; i<A.getLength(); i++)
		 {
		 cout << i << " (" << GetPosForIndexAtLevel(i, level).x << " " << GetPosForIndexAtLevel(i, level).y << ") " ;
		 grid[i].print();
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
			Coarsen(graph2, PQ, newIndex, unassigned, iNrOfCoarse, A);	
			if(bTiming) { SW.printTimeDiff();}
			
			if(bTiming) cout << endl;
		}
	}
#endif
	
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
	
#ifdef AMG_PRINT_COARSENING
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
	
#ifdef AMG_PRINT_P
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
	
#ifdef AMG_PRINT_R
	R.print();
#endif
	
	// create Galerkin product
	/////////////////////////////////////////	
	
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
	
#ifdef AMG_WRITE_MATRICES_PATH	
	if(A[0]->getLength() < AMG_WRITE_MATRICES_MAX)
	{
		cout << "write matrix A...";
		A[0]->writeToFile((string(AMG_WRITE_MATRICES_PATH) + "A" + nrstring(0) + ".mat").c_str());
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
	theta = 0.5;
	
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
	
	/*
	 cout << setw(level) << "AMG MGCycle, level " << level;
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
template<typename Matrix_type, typename Vector_type>
void amg<Matrix_type, Vector_type>::amgTestLevel(Vector_type &x, const Vector_type &b, int level)
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
void amg<Matrix_type, Vector_type>::amgTest(const Matrix_type& A_, Vector_type &vx, const Vector_type &b)
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
