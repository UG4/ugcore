/**
 * \file amg_debug.h
 *
 * \author Martin Rupp
 *
 * \date 16.06.10
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */



#ifndef __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_DEBUG_H__
#define __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_DEBUG_H__

#include "amg_nodeinfo.h"
#include "postscript.h"

namespace ug {
	
#if 0
//
// writeMatrices:
//----------------
//! writes A to pathAndName + "A" + level + ".mat for all levels, also R
//! @param pathAndName	path with name of matrix. for example "/Users/username/matrices/mat"
template<typename Matrix_type, typename Vector_type>
void amg<Matrix_type, Vector_type>::writeMatrices(const char *pathAndName)
{
	// only for small matrices
	if(A[0]->row_size() > 100*100*100*100)
		return;
	
	cout << "writing matrices "; cout.flush();
	string str(pathAndName);
	for(int i=0; i<used_levels-1; i++)
	{
		A[i]->writeToFile((str + "A" + ToString(i) + ".mat").c_str()); cout << "."; cout.flush();
		P[i].writeToFile((str + "P" + ToString(i) + ".mat").c_str()); cout << "."; cout.flush();
		R[i].writeToFile((str + "R" + ToString(i) + ".mat").c_str()); cout << "."; cout.flush();
	}
	if(used_levels > 0)
		A[used_levels-1]->writeToFile((str + "A" + ToString(used_levels-1) + ".mat").c_str());
	cout << " finished."; cout.flush();
}



// printCoarsening:
//----------------
//! Debug output. Writes position of Coarse nodes in coarse<level>.dat, and fine in fine<level>.dat for display in gnuplot
//! @param level	level which is to be printed
template<typename Matrix_type, typename Vector_type>
void amg<Matrix_type, Vector_type>::printCoarsening(int level, amg_nodeinfo *grid)
{  
	fstream fcoarse((string("/Users/mrupp/matrices/coarse") + ToString(level) + ".dat").c_str(), ios::out);
	fstream ffine  ((string("/Users/mrupp/matrices/fine") + ToString(level) + ".dat").c_str(), ios::out);
	int n = A[level]->row_size();
	
	for(int i=0; i < n; i++)
	{
		postype pos = GetPosForIndexAtLevel(i, level);
		if(grid[i].isCoarse())
			fcoarse << pos.x << " " << pos.y << " " << endl;
		else
			ffine << pos.x << " " << pos.y << " " << endl;
	}
  	
	/////////////
	
	fstream file((string("/Users/mrupp/matrices/coarsening") + ToString(level) + ".mat").c_str(), ios::out);
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
// amg<Matrix_type, Vector_type>::onlyOneLevel
//----------------
//! for testing. creates only one AMG level without direct solvers
//! @param A	matrix A.
template<typename Matrix_type, typename Vector_type>
bool amg<Matrix_type, Vector_type>::onlyOneLevel(const Matrix_type& A_)
{
	used_levels = 2;
	const Matrix_type *pA = &A_;
	A[0] = const_cast<Matrix_type*> (pA);
	
#ifdef AMG_WRITE_MATRICES_PATH	
	if(A[0]->row_size() < AMG_WRITE_MATRICES_MAX)
	{
		cout << "write matrix A...";
		A[0]->writeToFile((string(AMG_WRITE_MATRICES_PATH) + "A0.mat").c_str());
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
	
	cout << "[" << level << "]" << endl;
	
	
	stopwatch SW;
	SW.start();
	double pre1 = norm(b - Ah*x);
	cout << "calc norm "; SW.printTimeDiff(); cout << endl;
	
	SW.start();
	for(int i=0; i < nu1; i++)
		smoother[level].iterate(x, b);
	cout << nu1 << " presmoothing steps "; SW.printTimeDiff();
	
	Vector_type &r = *vec3[level];
	
	SW.start();	
	r = b - Ah*x;
	cout << "calculation of defect "; SW.printTimeDiff();
	
	double presmoothreduction = norm(r)/pre1;
	//cout << setw(2*level) << "presmooth reduction: " << norm(r)/pre1 << endl;
	
	Vector_type &rH = *vec1[level+1];
	Vector_type &eH = *vec2[level+1];
	
	SW.start();	
	rH = R[level]*r; 
	cout << "restriction "; SW.printTimeDiff();
	
	eH = 0.0;
	
	SW.start();	
	if(level+1 == used_levels-1)
		MGCycle(eH, rH, level+1);
	else
		for(int i=0; i<gamma; i++)
			MGCycle(eH, rH, level+1);
	cout << (level+1 == used_levels-1 ? 1 : gamma) << " MGCycles "; SW.printTimeDiff();
	
	x += P[level]*eH; //interpolate(eH, level);	
	double pre2 = norm(b - Ah*x);
	
	double res;
	
	SW.start();
	for(int i=0; i < nu2; i++)
		res = smoother[level].iterate(x, b);
	cout << nu2 << " postsmoothing steps "; SW.printTimeDiff();
	
	double post = norm(b - Ah*x);
	
	//cout << setw(2*level) << "postsmooth reduction: " << post/pre2 << endl;
	//cout << setw(2*level) << "level reduction: " << post/pre1 << endl;
	
	cout << endl;
	
	cout << "prered: " << presmoothreduction << " postred: " << post/pre2
	<< " levelred: " << post/pre1;
	
	if(post < 1e-10 || post/pre1 < 0.3)
		cout << " -> level seems to be OK." << endl;
	else
		cout << " -> LEVEL BROKEN!" << endl;
	
	
	cout << "norm(b-Ah*x) = " << post << endl;
	
	
	/*	vector<sortStruct<double> > inner;
	 vector<sortStruct<double> > boundary;
	 
	 for(int i=0; i<r.size(); i++)
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
	 */	
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
#endif // #if 0


// WriteToFile
//--------------------------------------------------
//! writes to a file in somewhat SparseMatrix-market format (for connection viewer)
template<typename T>
void AMGWriteToFile(const SparseMatrix<T> &A, int fromlevel, int tolevel, const char *filename, const cAMG_helper &h)
{
	fstream file(filename, ios::out);
	file << 1 << endl; // connection viewer version

	h.writePosToStream(file);
	file << 1 << endl;
	for(size_t i=0; i < A.num_rows(); i++)
	{
		for(typename SparseMatrix<T>::cRowIterator conn = A.beginRow(i); !conn.isEnd(); ++conn)
			if((*conn).dValue != 0.0)
				file << h.GetOriginalIndex(tolevel, i) << " " << h.GetOriginalIndex(fromlevel, (*conn).iIndex) << " " << ((*conn).dValue) << endl;
	}
}

// writeToFile
//--------------------------------------------------
//! writes to a file in somewhat SparseMatrix-market format (for connection viewer)
template<typename T>
void AMGWriteToFilePS(const SparseMatrix<T> &A, int fromlevel, int tolevel, const char *filename, const cAMG_helper &h)
{
	postscript ps;
	ps.create(filename);

	for(size_t i=0; i < A.num_rows(); i++)
	{
		int from = h.GetOriginalIndex(tolevel, i);
		ps.move_to(h.positions[from].x, h.positions[from].y);
		ps.print_text( string("0") + ToString(i));

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

}



#endif // __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_DEBUG_H__
