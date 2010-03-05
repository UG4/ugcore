/*
 *  experimental.h
 *  flexamg
 *
 *  Created by Martin Rupp on 18.02.10.
 *  Copyright 2010 G-CSC, University of Frankfurt. All rights reserved.
 *
 * unused code
 */


/////////////// FLENS /////////////////////////////////

#define VECLIB
#include <flens/flens.h>

using namespace flens;
using namespace std;


template<typename T>
int inverse(T &A)
{
	DenseVector<Array<int> > P(A.numRows());
	int info = trf(A, P);
	if(info) return info;
	info = tri(A, P);
	return info;
}
GeMatrix<FullStorage<double, ColMajor> > fSM( SM.getRows(), SM.getCols(), 0, 0);
for(int i=0; i<SM.getRows(); i++)
for(int j=0; j<SM.getCols(); j++)
{
	fSM(i,j) = SM(i,j);
}
cout << "_____________________________________________________________" << endl << endl;

cout << fSM;
GeMatrix<FullStorage<double, ColMajor> > fSM2(SM.getRows(), SM.getCols(), 0, 0); fSM2 = fSM;
inverse(fSM);

cout << fSM;





/*#define VECLIB
 #include <flens/flens.h>
 
 using namespace flens;
 using namespace std;
 
 
 template<typename T>
 int inverse(T &A)
 {
 DenseVector<Array<int> > P(A.numRows());
 int info = trf(A, P);
 if(info) return info;
 info = tri(A, P);
 return info;
 }*/



////////////////////////////// F-F Coarsening ///////////////////////////////////

template<typename Matrix_type, typename Vector_type>
int amg<Matrix_type, Vector_type>::Coarsen(cgraph &graph, maxheap<nodeinfo> &PQ, int *newIndex, int unassigned, int &iNrOfCoarse, const Matrix_type &A)
{
	Coarsen(graph, PQ, newIndex, unassigned, iNrOfCoarse, A);
	
#ifdef AMG_WRITE_COARSENING
	string s = (string(AMG_WRITE_MATRICES_PATH) + "A" + nrstring(A.tolevel) + ".mat");
	fstream fstr(s.c_str(), ios::out|ios::app);
#endif
	
	cout << iNrOfCoarse << " Coarse Nodes. ";
	
	// second pass
	//----------------
	// seems to work but doesnt help with convergence rates on complicated geometries???	
	cgraph TG;
	TG.createAsTransposeOf(graph);
	
	int nrOfFFCoarseNodes=0;
	int nrFF=0;
	
	vector<bool> marks(graph.getLength());
	
#if 0
	// somewhat more intelligent F-F prevention
	PQ.reset();	
	unassigned =0;
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
			if(grid[it()].isCoarse() || A.isUnconnected(i))
				continue;
			
			cgraph::cRowIterator it2(TG, it());
			for(; !it2.isEnd(); ++it2)
			{
				if(grid[it2()].isCoarse() && marks[it2()])
					break;
			}	
			
			if(it2.isEnd())
			{
#ifdef AMG_WRITE_COARSENING
				fstr << GetOriginalIndex(A.tolevel, i) << endl << GetOriginalIndex(A.tolevel, it()) << endl;
#endif
				// add both nodes again to coarsening process
				// rating = nr of F-F nodes involved.
				if(grid[i].rating < 0)
				{
					unassigned++;
					grid[i].rating = 1;
					PQ.insertItem(i);
				}
				else
				{
					grid[i].rating++;
					PQ.upheap(i);
				}
				
				
				int n = it();
				if(grid[n].rating < 0)
				{
					unassigned++;
					grid[n].rating = 1;
					PQ.insertItem(n);
				}
				else
				{
					grid[n].rating++;
					PQ.upheap(n);
				}
				nrFF++;
				
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
	
	if(unassigned == 0)
		return iNrOfCoarse;
	cout << unassigned << " nodes in F-F connections, coarsening... ".
	
	// Coarsen again with new ratings of F-F connections
#ifdef AMG_WRITE_COARSENING
	fstr.close();
#endif
	Coarsen(graph, PQ, newIndex, unassigned, iNrOfCoarse, A);
	cout << "now " << iNrOfCoarse << " coarse nodes. " << endl;
	
#ifdef AMG_WRITE_COARSENING
	fstr.open(s.c_str(), ios::out|ios::app);
#endif
#endif
	// third pass
	//----------------
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
			if(grid[it()].isCoarse() || A.isUnconnected(i))
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
				fstr << GetOriginalIndex(A.tolevel, i) << endl << GetOriginalIndex(A.tolevel, it()) << endl;
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
	return 0;
}