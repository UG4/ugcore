/*
 *  amg.hpp
 *  flexamg
 *
 *  Created by Martin Rupp on 03.12.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 *  implementation file for amg
 */
#pragma once
//#include "sparsematrix_util.h"

#include "amg_nodeinfo.h"
#include "stopwatch.h"


using namespace std;
namespace ug{


string nrstring(double d)
{
	char s[255];
	sprintf(s, "%g", d);
	return string(s);
}

string nrstring(int i)
{
	char s[255];
	sprintf(s, "%d", i);
	return string(s);
}

string nrstring(size_t i)
{
	char s[255];
	sprintf(s, "%u", (unsigned int) i);
	return string(s);
}

//#define GRAPH_WITH_LOCAL_INVERSE




#define AMG_WRITE_MATRICES_PATH "/Users/mrupp/matrices/AMG_"
#define AMG_WRITE_MATRICES_MAX (200*200)
#if 0


#define LATE_COARSE_SOLVER // do coarsening down to 10 nodes.
//#define AMG_WRITE_GRAPH


//#define AMG_PRINT_INDIRECT

#define AMG_PRINT_GRAPH

#define AMG_WRITE_COARSENING
#define AMG_WRITE_GRAPH


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


// WriteToFile
//--------------------------------------------------
//! writes to a file in somewhat SparseMatrix-market format (for connection viewer)
template<typename T>
void WriteToFile(const SparseMatrix<T> &A, int fromlevel, int tolevel, const char *filename, const cAMG_helper &h)
{
	fstream file(filename, ios::out);
	file << 1 << endl; // connection viewer version
	file << 2 << endl; // dimensions

	h.writePosToStream(file);
	file << 1 << endl;
	for(size_t i=0; i < A.num_rows(); i++)
	{
		for(typename SparseMatrix<T>::cRowIterator conn = A.beginRow(i); !conn.isEnd(); ++conn)
			if((*conn).dValue != 0.0)
				file << h.GetOriginalIndex(tolevel, i) << " " << h.GetOriginalIndex(fromlevel, (*conn).iIndex) << " " << ((*conn).dValue) << endl;
	}
}

class postscript
{
public:
	postscript()
	{
		bounds_left = 1e18;
		bounds_right = -1e18;
		bounds_top = 1e18;
		bounds_bottom = -1e18;
		mindist=1e18;
	}

	~postscript()
	{
		double scale = 1;
		scale = max( 400/(bounds_right-bounds_left), 400/(bounds_top-bounds_bottom) );
		file <<		"%%BoundingBox: " << bounds_left*scale*1.05 << " " << bounds_top*scale*1.05 << " " << bounds_right*scale*1.05 << " " << bounds_bottom*scale*1.05 << "\n"
					"%%Pages: 1\n"
					"%%DocumentsFonts: Monaco\n"
					"%%Copyright 2010 G-CSC - All Rights Reserved Worldwide\n"
					"%%EndComments\n\n";

		file <<		"1 setlinejoin\n"
					"1 setlinecap\n"
					"/Monaco findfont 10 scalefont setfont\n\n";

		file << 	"/M {moveto} def\n"
					"/S {lineto stroke} def\n"
					"/L {lineto} def\n"
					"/C {closepath fill} def\n"
					"/N {newpath} def\n"
					"/R {setrgbcolor} def\n"
					"/W {setlinewidth} def\n\n";

		file <<		"%%Endprolog\n"
					"%\n"
					"%%Page: 1 1\n"
					"%\n\n";

		mindist = sqrt(mindist);
		file << 	mindist * 0.01 << " W\n"
					"/Monaco findfont " << mindist * 0.2 << " scalefont setfont\n0 0 0 R\n";


		file << 	scale << " " << scale << " scale\n";
		file << out.str();
		file << "showpage\n\n%%Trailer";
	}

	bool create(const char *filename)
	{
		file.open((string(filename) + ".ps").c_str(), ios::out);
		file << 	"%!PS-Adobe-2.0 EPSF-1.2\n"
					"%%Title: " << filename << "\n"
					"%%Creator: ug postscript output\n"
					"%%CreationDate:\n";
		return true;
	}

	void setcolor(double r, double g, double b)
	{
		out << r << " " << g << " " << b << " R\n";
	}
	void move_to(double x, double y)
	{
		extend_bounds(x, y);
		out << "N " << x << " " << y << " M\n";
		last_movetox = x;
		last_movetoy = y;
	}

	void line_to(double x, double y)
	{
		extend_bounds(x, y);
		out << x << " " << y << " S\n";
		double d = ((x-last_movetox)*(x-last_movetox) + (y-last_movetoy)*(y-last_movetoy));
		if(d < mindist) mindist = d;
	}

	void extend_bounds(double x, double y)
	{
		if(x < bounds_left) bounds_left = x;
		if(x > bounds_right) bounds_right = x;
		if(y < bounds_top) bounds_top = y;
		if(y > bounds_bottom) bounds_bottom = y;
	}

	void line(double x1, double y1, double x2, double y2)
	{
		move_to(x1, y1);
		line_to(x2, y2);
	}

	void set_line_width(double width)
	{
		out << width << " W\n";
	}

	void print_text(const char *text)
	{
		out << "(" << text << ") show\n";
	}


private:
	double mindist;
	double last_movetox;
	double last_movetoy;
	double bounds_left, bounds_right, bounds_top, bounds_bottom;
	std::ostringstream out;

	fstream file;
};
// writeToFile
//--------------------------------------------------
//! writes to a file in somewhat SparseMatrix-market format (for connection viewer)
template<typename T>
void writeToFilePS(const SparseMatrix<T> &A, int fromlevel, int tolevel, const char *filename, const cAMG_helper &h)
{
	postscript ps;
	ps.create(filename);

	for(size_t i=0; i < A.num_rows(); i++)
	{
		int from = h.GetOriginalIndex(tolevel, i);
		ps.move_to(h.positions[from].x, h.positions[from].y);
		ps.print_text( (string("0") + nrstring(i)).c_str() );

		for(typename SparseMatrix<T>::cRowIterator conn = A.beginRow(i); !conn.isEnd(); ++conn)
		{
			if((*conn).dValue != 0.0)
			{
				if((*conn).iIndex != i)
				{
					int to = h.GetOriginalIndex(fromlevel, (*conn).iIndex);
					ps.move_to(h.positions[from].x, h.positions[from].y);
					ps.line_to(h.positions[to].x, h.positions[to].y);

				}
			}
		}
	}

	cout << endl;
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
void amg<Matrix_type, Vector_type>::CreateGraph(const Matrix_type &A, cgraph &graph, maxheap<amg_nodeinfo> &PQ, int &unassigned, amg_nodeinfo *grid)
{
	unassigned = 0;
	for(size_t i=0; i< A.num_rows(); i++)
	{
		graph.init(i);
		if(A[i].is_unconnected())
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
				graph.setConnection(i, (*conn).iIndex);
		}

		assert(graph.getNrOfConnections(i) > 0);
	}

#ifdef AMG_PRINT_GRAPH
	graph.print();
#endif

#ifdef AMG_WRITE_GRAPH
	graph.writeToFile((string(AMG_WRITE_MATRICES_PATH) + "G" + nrstring(A.tolevel) + ".mat").c_str(), amghelper, A.tolevel);
#endif

	// we need the transpose, since when we set a node coarse, we want
	// all nodes to be fine which can be interpolated by this coarse node
	// graph is afterwards made up of connections from a node i to j if
	// j has a strong connection to i
	graph.transpose();

#ifdef AMG_WRITE_GRAPH
	graph.writeToFile((string(AMG_WRITE_MATRICES_PATH) + "GT" + nrstring(A.tolevel) + ".mat").c_str(), A.tolevel);
#endif

	for(size_t i=0; i < A.num_rows(); i++)
	{
		if(A.is_unconnected(i))
			grid[i].setFineDirect();
		/*if(graph.getNrOfConnections(i) == 0)
			grid[i].setCoarse();*/
		else
		{
			//UG_ASSERT(graph.iNrOfConnections[i] > 0, "node " << i << " has " << graph.iNrOfConnections[i] << " connections?");
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
void amg<Matrix_type, Vector_type>::CreateGraph2(cgraph &graph, cgraph &graph2, maxheap<amg_nodeinfo> &PQ, int &unassigned, int &iNrOfCoarse, int *posInConnections, int *newIndex, amg_nodeinfo *grid)
{
	vector<int> connection(255);
	vector<int> nrOfPaths(255);

	amg_nodeinfo *grid2 = new amg_nodeinfo[graph.size()];

	iNrOfCoarse = 0;
	PQ.create(graph.size(), grid2);
	//graph.print();
	unassigned=0;
	for(size_t i=0; i < graph.size(); i++)
	{
		graph2.init(i);
		if(grid[i].isFineDirect())
		{
			grid2[i].setFineDirect();
			continue;
		}

		connection.clear();
		nrOfPaths.clear();
		// first calculate all nodes reachable with paths of length 2

		// ! i is coarse -> has only fine neighbors
		for(cgraph::cRowIterator conn(graph, i); !conn.isEnd(); ++conn)
		{
			size_t indexN = conn();
			for(cgraph::cRowIterator connN (graph, indexN); !connN.isEnd(); ++connN)
			{
				size_t indexNN = connN();

				if(indexNN == i || grid[indexNN].isFineDirect())
					continue;
				int pos = posInConnections[indexNN];
				if(pos == -1)
				{
					// never reached node indexNN from i, init.
					pos = posInConnections[indexNN]= connection.size();
					connection.push_back(indexNN);
					nrOfPaths.push_back(1);
				}
				else
					nrOfPaths[pos]++;
			}
		}

		// then sort out those which were reached #aggressiveCoarseningNrOfPaths (2 or 1) times
		grid2[i].rating = 0;
		for(size_t j=0; j<connection.size(); j++)
		{
			if(nrOfPaths[j] >= aggressiveCoarseningNrOfPaths)
			{
				// add connection i -> node
				graph2.setConnection(i, connection[j]);
				// increase rating of i
				grid2[i].rating ++;
			}
			// reset posInConnections for further use
			posInConnections[connection[j]] = -1;
		}

		// add node with rating > 0 to priority queue
		if(grid2[i].rating > 0)
		{
			PQ.insertItem(i);
			unassigned++;
		}
		else
		{
			grid2[i].setCoarse();
			newIndex[i] = iNrOfCoarse++;
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
int amg<Matrix_type, Vector_type>::Coarsen(cgraph &graph, maxheap<amg_nodeinfo> &PQ, int *newIndex, int unassigned, int &iNrOfCoarse, const Matrix_type &A, amg_nodeinfo *grid)
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
		UG_ASSERT(!grid[best].isAssigned(), "node " << best << " is already assigned??? (rating = " << grid[best].rating << ", unassigend = " << unassigned << ")");


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

	UG_ASSERT(iNrOfCoarse > 0, "no coarse nodes???");

	cout << iNrOfCoarse << " Coarse Nodes. ";

	// second pass
	//----------------
	// seems to work but doesnt help with convergence rates on complicated geometries???
	cgraph TG;
	TG.createAsTransposeOf(graph);

	int nrOfFFCoarseNodes=0;
	vector<bool> marks(graph.size());

	// seems to work but doesnt help with convergence rates on complicated geometries???

#ifdef AMG_WRITE_COARSENING
	fstr << "c" << endl;
#endif

	for(size_t i=0; i< graph.size(); i++)
	{
		if(grid[i].isCoarse() || A.is_unconnected(i))
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
			if(grid[it()].isCoarse() || A.is_unconnected(it()))
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
//! uses amg_nodeinfo *grid.
//!	@param P				Matrix P: here goes the calculated prolongation
//! @param	A				Matrix A: matrix for which to calculate prolongation on next level
//! @param	newIndex		newIndex of coarse Node i on next coarser level
//! @param	iNrOfCoarse		nr of coarse nodes on this level
template<typename Matrix_type, typename Vector_type>
void amg<Matrix_type, Vector_type>::CreateProlongation(SparseMatrix<double> &P, const Matrix_type &A, int *newIndex, int iNrOfCoarse, int &unassigned, amg_nodeinfo *grid)
{
	P.create(A.num_rows(), iNrOfCoarse);
#ifdef FLEXAMG
	P.fromlevel = A.fromlevel+1;
	P.tolevel = A.tolevel;
	P.name = "AMG:P";

#endif
	vector<SparseMatrix<double>::connection> con(255);
	SparseMatrix<double>::connection c;
	// DIRECT INTERPOLATION
	unassigned=0;

	for(size_t i=0; i < A.num_rows(); i++)
	{
		if(grid[i].isCoarse())
		{
			// a coarse node
			//P[i].initWithoutDiag();
			SparseMatrix<double>::connection con;
			con.iIndex = newIndex[i];  assert(newIndex[i] != -1);
			con.dValue = 1.0;
			P.set_matrix_row(i, &con, 1);
		}
		else if(A[i].is_unconnected())
		{
			//P[i].initWithoutDiag(); // boundary values need not to be prolongated
		}
		else if(grid[i].isFineDirect())
		{
			// a non-interpolated fine node. calculate interpolation weights

			// calc min off-diag-entry
			double dmax = 0, connValue, maxConnValue = 0;
			double sumNeighbors =0, sumInterpolatory=0;


			double diag = amg_diag_value(A.get_diag(i));

			for(typename Matrix_type::cRowIterator conn = A.beginRow(i); !conn.isEnd(); ++conn)
			{
				if((*conn).iIndex == i) continue; // skip diag
				connValue = amg_offdiag_value((*conn).dValue);

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
			//if(eps_truncation_of_interpolation > 0)  // Ruge/Stuebe A.7.2.4 truncation of interpolation
			//	barrier = min(theta*dmax, eps_truncation_of_interpolation*maxConnValue);
			//else
				barrier = theta*dmax;

			con.clear();

			for(typename Matrix_type::cRowIterator conn = A.beginRow(i); !conn.isEnd(); ++conn)
			{
				if((*conn).iIndex == i) continue; // skip diagonal
				if(!grid[(*conn).iIndex].isCoarse()) continue;

				connValue = amg_offdiag_value((*conn).dValue);
				if(connValue > barrier)
					continue;
				c.iIndex = newIndex[(*conn).iIndex];   assert(c.iIndex >= 0);
				c.dValue = connValue;

				con.push_back(c);
				sumInterpolatory += connValue;
			}

			// connections hinzufügen
			if(con.size() > 0)
			{
				double alpha = - (sumNeighbors / sumInterpolatory) / diag;
				for(size_t j=0; j<con.size(); j++)
					con[j].dValue *= alpha;

				//UG_ASSERT(con.size() > 0, "0 connections in point i = " << i << " ?");
				P.set_matrix_row(i, &con[0], con.size());
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
			UG_ASSERT(aggressiveCoarsening != 0, "no aggressive Coarsening but node " << i << " is fine and indirect??");
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
//! uses amg_nodeinfo *grid.
//!	@param P				Matrix P: here goes the calculated prolongation
//! @param	A				Matrix A: matrix for which to calculate prolongation on next level
//! @param	newIndex		newIndex of coarse Node i on next coarser level
//! @param	iNrOfCoarse		nr of coarse nodes on this level
//! @param posInConnections		array of size A.num_rows() for speedup of neighbor-neighbor-calculation inited with -1.
template<typename Matrix_type, typename Vector_type>
void amg<Matrix_type, Vector_type>::CreateIndirectProlongation(SparseMatrix<double> &P, const Matrix_type &A, int *newIndex, int *posInConnections, int unassigned, amg_nodeinfo *grid)
{
	UG_ASSERT(aggressiveCoarsening, "indirect interpolation only for aggressive coarsening");
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
		for(size_t i=0; i<A.num_rows() && unassigned > 0; i++)
		{
			if(!grid[i].isUnassignedFineIndirect() || A[i].is_unconnected())
				continue;

			double diag = amg_diag_value(A.get_diag(i));
			// calculate min offdiag-entry
			double dmax = 0;

			for(typename Matrix_type::cRowIterator conn = A.beginRow(i); !conn.isEnd(); ++conn)
			{
				if((*conn).iIndex == i) continue; // skip diagonal
				double connValue = amg_offdiag_value((*conn).dValue);
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

			for(typename Matrix_type::cRowIterator conn = A.beginRow(i); !conn.isEnd(); ++conn)
			{
				size_t indexN = (*conn).iIndex;
				if(indexN == i) continue; // skip diagonal

				// we dont want fine nodes which were indirectly interpolated in THIS pass
				if(grid[indexN].isFineIndirectLevel(pass))
					continue;
				// all interpolate neighbors are now from pass (pass-1) (otherwise makes no sense)

				double connValue = amg_offdiag_value((*conn).dValue);
				sumNeighbors += connValue;
				if(connValue > theta * dmax)
					continue;

				UG_ASSERT(!grid[indexN].isCoarse(), "Node " << i << " indirect, but neighbor " <<  indexN << " coarse?");

				// now we look from which nodes this fine node is interpolated from
				typename SparseMatrix<double>::rowIterator conn2 = P.beginRow(indexN); // !!! P
				for(; !conn2.isEnd(); ++conn2)
				{
					size_t indexNN = (*conn2).iIndex;
					int pos = posInConnections[indexNN];

					if(pos == -1)
					{
						pos = posInConnections[indexNN] = con2.size();
						c.iIndex = indexNN; assert(c.iIndex >= 0);

						AssignMult(c.dValue, connValue, (*conn2).dValue);
						con2.push_back(c);
						//nrOfPaths.push_back(1);
					}
					else
					{
						AddMult(con2[pos].dValue, connValue, (*conn2).dValue);
						//nrOfPaths[pos]++;
					}
				}
			}

			for(size_t j=0; j<con2.size(); j++)
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

			UG_ASSERT(sumInterpolatory != 0.0, " numerical unstable?");
			double alpha =  /*1/sumInterpolatory; */ - (sumNeighbors / sumInterpolatory)/diag;
			for(size_t j=0; j<con.size(); j++)
			{
				//cout << con[j].dValue << " - N:" << sumNeighbors << " I: " << sumInterpolatory << " alpha: " << alpha << ". " << con[j].dValue*alpha << " : " << A.get_diag(i) << endl;
				con[j].dValue *= alpha;
			}

			// connections hinzufügen
			P.set_matrix_row(i, &con[0], con.size());

		}

#ifndef NDEBUG
		if(unassigned == oldUnassigned)
		{
			cout << endl << "unassigned nodes left: " << endl;
			for(size_t i=0; i<A.num_rows(); i++)
			{
				if(grid[i].isUnassignedFineIndirect())
				   cout << i << " ";
			}
		}
		UG_ASSERT(unassigned != oldUnassigned, "Pass " << pass << ": Indirect Interpolation hangs at " << unassigned << " unassigned nodes.")
#endif

#ifdef AMG_PRINT_INDIRECT
		cout << "calculated, ";
#endif
		cout << unassigned << " left. ";
		pass++;
		oldUnassigned = unassigned;
		//break;
	}

	P.finalize();
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
	amg_nodeinfo *grid = new amg_nodeinfo[A.num_rows()];

	bool bTiming=true;
	cout << "Creating level " << level << ". (" << A.num_rows() << " nodes)" << endl;
	stopwatch SW;
	stopwatch SWwhole; SWwhole.start();

	// amg_nodeinfo: infos zu den einzelnen knoten von A, verwaltung der ratings etc
	// entries: speicherung der strong connections, nur wohin
	// entriesvalues: speicherung der WERTE der connections.
	// (Werte werden bei Konstruktion des Grobgitters nicht benötigt)

	// todo: check for isolated condition

	maxheap<amg_nodeinfo> PQ(A.num_rows(), grid);

	std::vector<bool>	coarse(A.num_rows());
	int unassigned = A.num_rows();

	int *newIndex = new int[A.num_rows()];
	memset(newIndex, -1, sizeof(int)*A.num_rows());

	int *posInConnections = new int[A.num_rows()];
	memset(posInConnections, -1, sizeof(int)*A.num_rows());
	cout << "A.totalNrOfConnections = " << A.total_num_connections() << endl;
	cgraph graph(A.num_rows());


	// build graph
	/////////////////////////////////////////

	cout << "building graph... "; cout.flush();
	if(bTiming) SW.start();
	// CreateGraph(C, graph, PQ, unassigned);
	CreateGraph(A, graph, PQ, unassigned, grid);
	if(bTiming) SW.printTimeDiff();

#ifdef AMG_PRINT_COARSEN_RATINGS
	 for(size_t i=0; i<A.num_rows(); i++)
		 cout << i << " [" << GetOriginalIndex(level, i) << "] " << grid[i] << endl;
#endif
	//PQ.print();

	// Coarsen
	/////////////////////////////////////////
	cout << "coarsening... "; cout.flush();

	if(bTiming) SW.start();
	int iNrOfCoarse = 0;
	Coarsen(graph, PQ, newIndex, unassigned, iNrOfCoarse, A, grid);

	if(bTiming) { SW.printTimeDiff();}

	// agressive Coarsening
	/////////////////////////////////////////
#ifndef	GRAPH_WITH_LOCAL_INVERSE
	if(aggressiveCoarsening && level == 0)
	{
		// build graph 2
		//------------------

		cgraph graph2(A.num_rows());

		cout << "building graph2... "; cout.flush();
		if(bTiming) SW.start();

		memset(newIndex, -1, sizeof(int)*A.num_rows());

		CreateGraph2(graph, graph2, PQ, unassigned, iNrOfCoarse, posInConnections, newIndex, grid);
		if(bTiming) SW.printTimeDiff();

		/*for(int i=0; i<A.num_rows(); i++)
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
			Coarsen(graph2, PQ, newIndex, unassigned, iNrOfCoarse, A, grid);
			if(bTiming) { SW.printTimeDiff();}

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
	for(size_t i=0; i<A.num_rows(); i++)
		if(grid[i].isCoarse())
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
		if(grid[i].isCoarse())
			parentIndex[level+1][ newIndex[i] ] = i;
//#endif

	/*for(size_t i=0; i<A.num_rows(); i++)
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

#if FLEXAMG
	P.fromlevel = level+1;
	P.tolevel = level;
#endif

	if(bTiming) SW.start();
	CreateProlongation(P, A, newIndex, iNrOfCoarse, unassigned, grid);
	if(unassigned > 0)
		CreateIndirectProlongation(P, A, newIndex, posInConnections, unassigned, grid);

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
	CreateAsMultiplyOf(AH, R, A, P, posInConnections);
	//createGalerkinMatrix(AH, R, A, P, posInConnections);

	/*	SparseMatrix<double> RA;
	 RA.createAsMultiplyOf(R, A);
	 AH.createAsMultiplyOf(RA, P); */

	/*	SparseMatrix<double> AP;
	 AP.createAsMultiplyOf(A, P);
	 AH.createAsMultiplyOf(R, AP);*/


	if(bTiming) SW.printTimeDiff();
	AH.name = "AMG:A";
#ifdef FLEXAMG
	AH.fromlevel = level+1;
	AH.tolevel = level+1;
#endif

	if(bTiming) SW.start();
	AH.finalize();
	if(bTiming) { cout << "finalizing..."; SW.printTimeDiff(); }

#ifdef AMG_PRINT_AH
	cout << "AH level " << level << endl;
	AH.print();
#endif

	// finish
	/////////////////////////////////////////

	//AH.print("AH");
	int nnz = AH.total_num_connections();
	cout << "AH: nnz: " << nnz << " Density: " << double(nnz)/(double(AH.num_rows())*double(AH.num_rows()))*100.0 << "% nnz/n: " << nnz/(double)AH.num_rows() << endl;
	cout << "Coarsening rate: " << (100.0*AH.num_rows())/(A.num_rows()) << "%" << endl;

	cout << " level "; SWwhole.printTimeDiff();  cout << endl;
	cout << endl;
	cout.flush();

//	cout << A << endl;
//	cout << P << endl;
//	cout << R << endl;


#ifdef AMG_WRITE_MATRICES_PATH
	if(this->A[0]->num_rows() < AMG_WRITE_MATRICES_MAX)
	{
		cout << "write matrices";
		writeToFile(P, level+1, level, (string(AMG_WRITE_MATRICES_PATH) + "P" + nrstring(level) + ".mat").c_str(), amghelper); cout << "."; cout.flush();
		writeToFile(R, level, level+1, (string(AMG_WRITE_MATRICES_PATH) + "R" + nrstring(level) + ".mat").c_str(), amghelper); cout << "."; cout.flush();
		writeToFile(AH, level+1, level+1, (string(AMG_WRITE_MATRICES_PATH) + "A" + nrstring(level+1) + ".mat").c_str(), amghelper); cout << "."; cout.flush();
		cout << " done." << endl;
	}
#endif

	//P.print();

	cout << AH << endl;
	cout << endl;//*/

	delete[] posInConnections;
	delete[] newIndex;
	delete[] grid;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// amg<Matrix_type, Vector_type>::init
//----------------
//! creates MG Hierachy for with Matrix_type A and temporary vectors for higher levels
//! @param A	matrix A.
template<typename Matrix_type, typename Vector_type>
bool amg<Matrix_type, Vector_type>::init(const Matrix_type& A_)
{
	amghelper.positions = positions;
	amghelper.size = A_.num_rows();
	amghelper.parentIndex = parentIndex;

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
		writeToFile(*A[0], 0, 0, (string(AMG_WRITE_MATRICES_PATH) + "A" + nrstring(0) + ".mat").c_str(), amghelper);
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

		vec3[i] = new Vector_type (A[i]->num_rows(), "AMG:tempvec 3");
		vec3[i]->level = i;
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

	cout << "Operator Complexity: " << nnzs/A[0]->total_num_connections() << " grid complexity: " << totallength/A[0]->num_rows() << endl << endl;

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

	// apply lmgc on coarser grid
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


#ifdef FLEXAMG
#include "amg_debug.h"
#endif
