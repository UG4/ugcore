/*
 *  SparseMatrix.hpp
 *  flexamg
 *
 *  Created by Martin Rupp on 04.11.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 */
#pragma once
#pragma mark creation etc

#include <fstream>
////////////////////////////////////////////////////////////////////////////////
static int GetOriginalIndex(int level, int i) { return i; }

namespace ug{
//!
//! constructor for empty SparseMatrix
template<typename T>
SparseMatrix<T>::SparseMatrix()
{
	name = "?";
	cols = rows = iTotalNrOfConnections = 0;
	iFragmentedMem = 0;
	cons = NULL;
	consmem = NULL; consmemsize = 0;
	iTotalNrOfConnections = 0;
	iNrOfConnections = NULL;
	fromlevel = tolevel = -1;
	bandwidth = 0;
	
	estimatedRowSize = 0;
	iMaxNrOfConnections = NULL;
	
	FORCE_CREATION { print(); p(); pr(0); }
}
	
template<typename T> bool
SparseMatrix<T>::destroy()
{
  if(cons)
	for(int i=0; i < rows; i++)
		safeSetConnections(i, NULL);

  if(iNrOfConnections) delete [] iNrOfConnections; iNrOfConnections = NULL;
	if(consmem)	delete [] consmem; consmem = NULL;
	if(cons) delete [] cons; cons = NULL;
	return true;
}
	
	
//!
//! destructor
template<typename T>
SparseMatrix<T>::~SparseMatrix()
{
	destroy();
}

// create
//---------------------
//! used to create the SparseMatrix
//! @param _rows nr of rows
//! @param _cols nr of cols
template<typename T>
bool SparseMatrix<T>::create(int _rows, int _cols)
{
	UG_ASSERT(rows == 0 && cols == 0, *this << " not empty.");
	
	rows = _rows;
	cols = _cols;
	
	cons = new connection*[rows+1];
	memset(cons, 0, sizeof(connection*)*(rows+1));

	iNrOfConnections = new int[rows];
	memset(iNrOfConnections, 0, sizeof(int)*rows);
	
	iMaxNrOfConnections = new int[rows];
	memset(iMaxNrOfConnections, 0, sizeof(int)*rows);
	
	iTotalNrOfConnections = 0;
	bandwidth = 0;
	return true;
}

// createAsTransposeOf
//-----------------------
//!
//! write in a empty SparseMatrix the transpose SparseMatrix of B.
template<typename T>
void SparseMatrix<T>::createAsTransposeOf(const SparseMatrix &B)
{
	ASSERT1(B.cols > 0 && B.rows > 0);
	create(B.cols, B.rows);
	fromlevel = B.tolevel;
	tolevel = B.fromlevel;
	
	// get length of each row
	for(int j=0; j < B.rows; j++)
		for(cRowIterator conn = B.beginRow(j); !conn.isEnd(); ++conn)
			if((*conn).dValue == 0) continue;
			else iNrOfConnections[(*conn).iIndex]++;
	
	int newTotal = 0;
	for(int i=0; i < rows; i++)
		newTotal += iNrOfConnections[i];
	
	
	
	int *nr = new int[rows];
	// init SparseMatrix data structure
	consmem = new connection[newTotal];
	consmemsize = newTotal;
	
	connection *p = consmem;
	for(int i=0; i < rows; i++)
	{
		nr[i] = 0;
		cons[i] = p;
		p += iNrOfConnections[i];
		iMaxNrOfConnections[i] = iNrOfConnections[i];
	}
	
	iTotalNrOfConnections = newTotal;
	iFragmentedMem = 0;
	
	bandwidth = 0;
	// write values
	for(int i=0; i < B.rows; i++)
		for(cRowIterator conn = B.beginRow(i); !conn.isEnd(); ++conn)
		{
			if((*conn).dValue == 0) continue;
			int ndx = (*conn).iIndex;
			UG_ASSERT(ndx >= 0 && ndx < rows, "connection " << (*conn) << " of " << B << ", row " << i << " out of range 0.." << rows);
			
			int k= nr[ndx];
			
			UG_ASSERT(k>=0 && k<iNrOfConnections[ndx], "k = " << k << ". precalculated nr of Connections " << iNrOfConnections[ndx] << " wrong?");
			//ASSERT(cons[ndx] + k < consmem+iMaxTotalNrOfConnections);
			//ASSERT(cons[ndx] + k >= consmem);
			
			cons[ndx][k].dValue = (*conn).dValue;
			cons[ndx][k].iIndex = i;
			if(bandwidth < i-ndx) bandwidth = i-ndx;
			nr[ndx]++;
		}
	
	delete[] nr;
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// createAsMultiplyOf:
//-------------------------
//! Calculates *this = A B C. posInConnections only needed for speedup (has to be -1 forall i).
//! other possibility: search in vector con for con[i].iIndex == indexTo about 3-4x slower (3d)
//! @param  A
//! @param  B
//! @param  C
//! @param posInConnections		array of size B.getCols() for speedup of neighbor-neighbor-calculation inited with -1.
template<typename T>
template<typename A_type, typename B_type, typename C_type>
void SparseMatrix<T>::createAsMultiplyOf(const A_type &A, const B_type &B, const C_type &C, int *posInConnections)
{
	UG_ASSERT(C.getRows() == B.getCols() && B.getRows() == A.getCols(), "-");
	//cout << endl << " Creating Galerkin Matrix_type..." << endl;
	
	
	// speedup with array posInConnections, needs n memory
	// posInConnections[i]: index in the connections for current row.
	// has to be -1 for all nodes
	
	bool bOwnMem = false;
	if(posInConnections == NULL)
	{		
		posInConnections = new int[C.getCols()];
		for(int i=0; i<C.getCols(); i++) posInConnections[i] = -1; // memset(posInConnections, -1, sizeof(int)*C.getCols());
		bOwnMem = true;
		UG_ASSERT(0, "");
	}
#ifdef DEBUG
	else
	{
		for(int i=0; i<C.getCols(); i++) 
			UG_ASSERT(posInConnections[i] == -1, "posInConnections[" << i << "] has to be -1, but is " << posInConnections[i] << ".");
	}
#endif
	
	vector<connection > con(255);
	
	//typename A_type::entry_type a;
	//typename Mult_Traits<typename A_type::entry_type, typename B_type::entry_type>::ReturnType ab;
	//typename C_type::entry_type cvalue;
	double a, ab, cvalue;
	
	
	connection c;
	
	create(A.getRows(), C.getCols());	
	
	for(int i=0; i < A.getRows(); i++)
	{
		// we want to have the diagonal first:
		posInConnections[i] = 0;
		con.clear();
		c.iIndex = i;
		c.dValue = 0.0;
		con.push_back(c);
		
		
		for(typename A_type::cRowIterator itA(A, i); !itA.isEnd(); ++itA)
		{			
			if((*itA).dValue == 0.0) continue;
			a = (*itA).dValue;
			
			for(typename B_type::cRowIterator itB(B, (*itA).iIndex); !itB.isEnd(); ++itB)
			{
				if((*itB).dValue == 0.0) continue;
				assign_mult(ab, a, (*itB).dValue);
				
				for(typename C_type::cRowIterator itC(C, (*itB).iIndex); !itC.isEnd(); ++itC)
				{
					cvalue = (*itC).dValue;
					if(cvalue == 0.0) continue;					
					int indexTo = (*itC).iIndex;					
										
					if(posInConnections[indexTo] == -1)
					{
						// we havent visited node <indexTo>
						// so we need to add a Connection to the row
						// save the index of the connection in the row
						posInConnections[indexTo] = con.size();
						c.iIndex = indexTo;
						assign_mult(c.dValue, ab, cvalue);
						con.push_back(c);				
					}
					else
					{
						// we have visited this node before,
						// so we know the index of the connection
						// -> add a*b*c
						//TODO 
						add_mult(con[posInConnections[indexTo]].dValue, ab, cvalue);
					}
					
				}
			}
		}
		
		// reset posInConnections to -1
		for(int j=0; j<con.size(); j++) posInConnections[con[j].iIndex] = -1;				
		// set Matrix_type Row in AH
		setMatrixRow(i, &con[0], con.size());
	}
	
	if(bOwnMem)
		delete[] posInConnections;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// createAsMultiplyOf:
//-------------------------
//! Calculates *this = A B. posInConnections only needed for speedup (has to be -1 forall i).
//! other possibility: search in vector con for con[i].iIndex == indexTo about 3-4x slower (3d)
//! @param  A
//! @param  B
//! @param posInConnections		array of size B.getCols() for speedup of neighbor-neighbor-calculation inited with -1.
template<typename T>
template<typename A_type, typename B_type>
void SparseMatrix<T>::createAsMultiplyOf(const A_type &A, const B_type &B, int *posInConnections)
{
	UG_ASSERT(B.getRows() == A.getCols(), "");
	//cout << endl << " Creating Galerkin Matrix_type..." << endl;
	
	
	// speedup with array posInConnections, needs n memory
	// posInConnections[i]: index in the connections for current row.
	// has to be -1 for all nodes
	
	bool bOwnMem = false;
	if(posInConnections == NULL)
	{		
		posInConnections = new int[B.getCols()];
		for(int i=0; i<B.getCols(); i++) posInConnections[i] = -1; // memset(posInConnections, -1, sizeof(int)*C.getCols());
		bOwnMem = true;
	}
#ifdef DEBUG
	else
	{
		for(int i=0; i<B.getCols(); i++) 
			UG_ASSERT(posInConnections[i] == -1, "posInConnections[" << i << "] has to be -1, but is " << posInConnections[i] << ".");
	}
#endif
	
	vector<connection > con(255);
	
	typename A_type::entry_type a;
	typename B_type::entry_type b;
	
	connection c;	
	create(A.getRows(), B.getCols());	
	
	for(int i=0; i < A.getRows(); i++)
	{
		// we want to have the diagonal first:
		posInConnections[i] = 0;
		con.clear();
		c.iIndex = i;
		c.dValue = 0.0;
		con.push_back(c);
		
		for(typename A_type::cRowIterator itA(A, i); !itA.isEnd(); ++itA)
		{			
			if((*itA).dValue == 0.0) continue;
			a = (*itA).dValue;
			
			for(typename B_type::cRowIterator itB(B, (*itA).iIndex); !itB.isEnd(); ++itB)
			{
				if((*itB).dValue == 0.0) continue;
				b = (*itB).dValue; 
				
				int indexTo = (*itB).iIndex;					
				if(posInConnections[indexTo] == -1)
				{
					// we havent visited node <indexTo>
					// so we need to add a Connection to the row
					// save the index of the connection in the row
					posInConnections[indexTo] = con.size();
					c.iIndex = indexTo;
					c.dValue = a * b;
					con.push_back(c);				
				}
				else
				{
					// we have visited this node before,
					// so we know the index of the connection
					// -> add a*b*c
					//TODO 
					con[posInConnections[indexTo]].dValue += a * b;
				}
			}
		}
		
		// reset posInConnections to -1
		for(int j=0; j<con.size(); j++) posInConnections[con[j].iIndex] = -1;				
		// set Matrix_type Row in AH
		setMatrixRow(i, &con[0], con.size());
	}
	
	if(bOwnMem)
		delete[] posInConnections;
}


/*void SparseMatrix<T>::recreateWithMaxNrOfConnections(int newMax) const
 {
 // create new cons Memory
 connection *consmemNew = new connection[newMax];
 
 // adjust pointers
 int diff = consmemNew - consmem;
 for(int i=0; i<rows; i++)
 {
 if(cons[i] != NULL) 
 cons[i] += diff;
 }
 
 // copy, delete old and swap
 memcpy(consmemNew, consmem, sizeof(connection)*iTotalNrOfConnections);
 delete[] consmem;
 consmem = consmemNew;
 iMaxTotalNrOfConnections = newMax;
 }*/



#pragma mark general functions
////////////////////////////////////////////////////////////////////////////////



// eliminateDirichletValues
//----------------------------
//!
//! eliminates Dirichlet Values by putting them on the rhs b and setting the row to i,i = 1.0;
template<typename T>
template<typename Vector_type>
void SparseMatrix<T>::eliminateDirichletValues(Vector_type &b)
{
	for(int i=0; i<rows; i++)
	{
		if(isUnconnected(i)) continue;		
		for(rowIterator conn = beginRow(i); !conn.isEnd(); ++conn)
		{
			int conindex = (*conn).iIndex;
			if(isUnconnected(conindex))
			{
				sub_mult(b[i], (*conn).dValue, b[conindex]);
				(*conn).dValue = 0;
			}
		}
	}
}

// setDirichletRow
//----------------------------
//!
//! sets the row to i,i = 1.0.
template<typename T>
void SparseMatrix<T>::setDirichletRow(int row)
{
	UG_ASSERT(row >= 0 && row < rows, *this << ": row " << row << " out of bounds.");
	if(iNrOfConnections[row] == 0)
	{
		connection *c = new connection[1];
		safeSetConnections(row, c);
	}
	cons[row][0].dValue = 1.0;
	cons[row][0].iIndex= row;

	iNrOfConnections[row] = 1;
	iMaxNrOfConnections[row] = 1;
}

//!
//! sets the # nrows in pRows to i,i = 1.0.
template<typename T>
void SparseMatrix<T>::setDirichletRows(int *pRows, int nrows)
{
	for(int i=0; i<nrows; i++)
		setDirichletRow(pRows[i]);
}

template<typename T>
bool SparseMatrix<T>::set_dirichlet_rows(const local_index_type &ind)
{
	for(std::size_t i=0; i<ind.size(); i++)
		setDirichletRow(ind[i][0]);
	return true;
}

// res = A*x
template<typename T>
template<typename Vector_type>
bool SparseMatrix<T>::apply(Vector_type &res, const Vector_type &x) const
{
	UG_ASSERT(cols == x.size(), "x: " << x << " has wrong length (should be " << cols << "). A: " << *this);
	UG_ASSERT(rows == res.size(), "res: " << x << " has wrong length (should be " << rows << "). A: " << *this);
	
	/*
	int n = min(rows, cols);
	for(int i=0; i < n; i++)
	{
		cRowIterator conn = beginRow(i);
		res[i] = 0.0;
		for(; !conn.isEnd(); ++conn)
			res[i] += (*conn).dValue * x[(*conn).iIndex];
	}
	
	for(int i=n; n < cols; i++)
	{
		cRowIterator conn = beginRow(i);
		++conn; // skip diag, since x.size <= n
		res[i] = 0.0;
		for(; !conn.isEnd(); ++conn)
			res[i] += (*conn).dValue * x[(*conn).iIndex];
	}*/

	for(int i=0; i < rows; i++)
	{
		cRowIterator conn = beginRow(i);
		res[i] = 0.0;

		for(; !conn.isEnd(); ++conn)
		{
			if((*conn).dValue != 0)
				res[i] += (*conn).dValue * x[(*conn).iIndex];
		}
	}

	return true;
}



// res = A.T() * x
template<typename T>
template<typename Vector_type>
bool SparseMatrix<T>::applyTransposed(Vector_type &res, const Vector_type &x) const
{
	UG_ASSERT(rows == x.size(), "x: " << x << " has wrong length (should be " << rows << "). A: " << *this);
	UG_ASSERT(cols == res.size(), "res: " << x << " has wrong length (should be " << cols << "). A: " << *this);
	
	res = 0.0;
	/*
	int n = min(rows, x.size());
	for(int i=0; i<n; i++)
	{
		cRowIterator conn = beginRow(i);
		for(; !conn.isEnd(); ++conn)
			res[(*conn).iIndex] += (*conn).dValue * x[i];
	}
	*/

	/*int n = min(cols, rows);
	for(int i=0; i<n; i++)
	{
		cRowIterator conn = beginRow(i);
		for(; !conn.isEnd(); ++conn)
			res[(*conn).iIndex] += (*conn).dValue * x[i];
	}

	for(int i=n; i<rows; i++)
	{
		cRowIterator conn = beginRow(i);
		++conn;
		for(; !conn.isEnd(); ++conn)
			res[(*conn).iIndex] += (*conn).dValue * x[i];
	}*/
	for(int i=0; i<rows; i++)
		{
			cRowIterator conn = beginRow(i);
			++conn;
			for(; !conn.isEnd(); ++conn)
			{
				if((*conn).dValue != 0.0)
					res[(*conn).iIndex] += (*conn).dValue * x[i];
			}
		}
	// ??? since x.size() <= n, x[i] doesnt exist for i >=n.
	return true;
}


// res = res - A*x
template<typename T>
template<typename Vector_type>
bool SparseMatrix<T>::matmul_minus(Vector_type &res, const Vector_type &x) const
{
	UG_ASSERT(cols == x.size(), "x: " << x << " has wrong length (should be " << cols << "). A: " << *this);
	UG_ASSERT(rows == res.size(), "res: " << x << " has wrong length (should be " << rows << "). A: " << *this);
	
	int n = min(rows, x.size());
	for(int i=0; i < n; i++)
	{
		cRowIterator conn = beginRow(i);
		for(; !conn.isEnd(); ++conn)
			res[i] -= (*conn).dValue * x[(*conn).iIndex];
	}
	
	for(int i=n; n < rows; i++)
	{
		cRowIterator conn = beginRow(i);
		++conn; // skip diag, since x.size <= n
		for(; !conn.isEnd(); ++conn)
			res[i] -= (*conn).dValue * x[(*conn).iIndex];
	}
	return true;
}


//!
//! @return a matrixrow object of row i
template<typename T>
const matrixrow<T> SparseMatrix<T>::getrow(int i) const
{	
	return matrixrow<entry_type> (*this, i);
}

//!
//! @return a matrixrow object of row i
template<typename T>
const matrixrow<T> SparseMatrix<T>::operator [] (int i) const
{		
	return matrixrow<entry_type> (*this, i);
}


// getDiag
//-------------
//! get Diagonal A_[i,i] of matrix
//! @param i
//! @return A_{i,i}
template<typename T>
inline const T &SparseMatrix<T>::getDiag(int i) const
{
	UG_ASSERT(cons[i][0].iIndex == i, *this << " first entry has to be diagonal");
	// evtl anders, da nicht jede Matrix diageintrag
	return cons[i][0].dValue;
}

template<typename T>
inline T &SparseMatrix<T>::getDiag(int i)
{
	UG_ASSERT(cons[i][0].iIndex == i, *this << " first entry has to be diagonal");
	return cons[i][0].dValue;
}

// isUnconnected
//-----------------
/*!
 @remark since first connection is always i,i, every row has at least 1 connection
 @return true if row i has no connection to indices other than i
 */
template<typename T>
inline bool SparseMatrix<T>::isUnconnected(int i) const
{
	UG_ASSERT(i < rows && i >= 0, *this << ": " << i << " out of bounds.");
	return iNrOfConnections[i] <= 1; 
}




// add Submatrix
//--------------------
//! function to add submatrices ( submatrix )
template<typename T>
template<typename M>
void SparseMatrix<T>::add(const M &mat, int *rows, int *cols)
{
	connection *c = new connection[mat.getCols()];
	int nc;
	for(int i=0; i < mat.getRows(); i++)
	{
		nc = 0;
		for(int j=0; j < mat.getCols(); j++)
		{
			if(mat(i,j) != 0.0)
			{
				c[nc].iIndex = cols[j];
				c[nc].dValue = mat(i, j);
				nc++;
			}
		}
		UG_ASSERT(nc < mat.getCols(), "???");

		if(nc > 0)
			addMatrixRow(rows[i], c, nc);
	}
	delete[] c;
}

// set Submatrix
//! 
//! function to add submatrices ( \sa submatrix )
template<typename T>
template<typename M>
void SparseMatrix<T>::set(const M &mat, int *rows, int *cols)
{
	connection *c = new connection[mat.getCols()];
	int nc;
	for(int i=0; i < mat.getRows(); i++)
	{
		nc = 0;
		for(int j=0; j < mat.getCols(); j++)
		{
			if(mat(i,j) != 0.0)
			{
				c[nc].iIndex = cols[j];
				c[nc].dValue = mat(i, j);
				nc++;
			}
		}
		UG_ASSERT(nc < mat.getCols(), "???");
		if(nc > 0)
			setMatrixRow(rows[i], c, nc);
	}
	delete[] c;
}

//!
template<typename T>
template<typename M>
void SparseMatrix<T>::get(M &mat, int *rows, int *cols) const
{
	vector<sortStruct<int> > sortedCols(mat.getCols());
	
	for(int i=0; i<mat.getCols(); i++)
	{
		sortedCols[i].index = i;
		sortedCols[i].sortValue = cols[i];
	}	
	sort(sortedCols.begin(), sortedCols.end());
	
	for(int i=0; i < mat.getRows(); i++)
	{
		int iindex_global = rows[i];
		int iindex_local = i;
	
		cRowIterator conn = beginRow(iindex_global);
		// diagonal
		mat(iindex_local, iindex_local) = (*conn).dValue;
		++conn;

		int j = 0;
		while(j < mat.getCols() && !conn.isEnd())
		{
			int jindex_global = sortedCols[j].sortValue;
			int jindex_local = sortedCols[j].index;

			if(iindex_global == jindex_global) { j++; continue; }
			int cindex_global = (*conn).iIndex;
			
			if(cindex_global < jindex_global) 
				++conn;
			else if(cindex_global > jindex_global) 
			{
				mat(iindex_local, jindex_local) = 0.0;
				++j;
			}
			else
			{
				mat(iindex_local, jindex_local) = (*conn).dValue;
				++conn; ++j;
			}
		}
	}
}

template<typename T>
template<typename M>
void SparseMatrix<T>::add(const M &mat, vector<int> &rows, vector<int> &cols)
{
	ASSERT1(mat.getCols() == cols.size() && mat.getRows() == rows.size());
	add(mat, &rows[0], &cols[0]);
}

template<typename T>
template<typename M>
void SparseMatrix<T>::set(const M &mat, vector<int> &rows, vector<int> &cols)
{
	ASSERT1(mat.getCols() == cols.size() && mat.getRows() == rows.size());
	set(mat, &rows[0], &cols[0]);
}

template<typename T>
template<typename M>
void SparseMatrix<T>::get(M &mat, vector<int> &rows, vector<int> &cols) const
{
	ASSERT1(mat.getCols() == cols.size() && mat.getRows() == rows.size());
	get(mat, &rows[0], &cols[0]);
}

//////////////////////////////////////////////////////////////////////////
template<typename T>
bool SparseMatrix<T>::get(local_matrix_type& mat, const local_index_type& I, const local_index_type& J) const
{
	for(uint i = 0; i < I.size(); ++i)
	{
		// do this sorted
		for(uint j = 0; j < J.size(); ++j)
		{
			mat(i,j) = 0.0;
			for(cRowIterator conn = beginRow(I[i][0]); !conn.isEnd(); ++conn)
			{
				if((*conn).iIndex == J[j][0])
				{
					mat(i,j) = (*conn).dValue;
					break;
				}
			}
		}
	}
	return true;
}

template<typename T>
bool SparseMatrix<T>::add(const local_matrix_type& mat, const local_index_type& I, const local_index_type& J)
{
	connection *c = new connection[J.size()];
	std::size_t nc;
	for(std::size_t i=0; i < I.size(); i++)
	{
		nc = 0;
		for(std::size_t j=0; j < J.size(); j++)
		{
			if(mat(i,j) != 0.0)
			{
				c[nc].iIndex = J[j][0];
				c[nc].dValue = mat(i, j);
				nc++;
			}
		}
		UG_ASSERT(nc <= J.size(), "???");
		if(nc > 0)
			addMatrixRow(I[i][0], c, nc);
	}

	delete[] c;
	return true;
}

template<typename T> bool SparseMatrix<T>::set(const local_matrix_type& mat, const local_index_type& I, const local_index_type& J)
{
	connection *c = new connection[J.size()];
	std::size_t nc;
	for(int i=0; i < I.size(); i++)
	{
		nc = 0;
		for(int j=0; j < J.size(); j++)
		{
			if(mat(i,j) != 0.0)
			{
				c[nc].iIndex = J[j][0];
				c[nc].dValue = mat(i, j);
				nc++;
			}
		}
		UG_ASSERT(nc <= J.size(), "???");
		if(nc > 0)
			setMatrixRow(I[i][0], c, nc);
	}
	delete[] c;
	return true;
}
//////////////////////////////////////////////////////////////////////////


template<typename T>
void SparseMatrix<T>::add(const entry_type &d, int row, int col)
{
	connection c;
	c.dValue = d;
	c.iIndex = col;
	addMatrixRow(row, &c, 1);
}

template<typename T>
void SparseMatrix<T>::set(const entry_type &d, int row, int col)
{
	connection c;
	c.dValue = d;
	c.iIndex = col;
	setMatrixRow(row, &c, 1);
}

template<typename T>
void SparseMatrix<T>::get(entry_type &d, int row, int col) const
{
	if(row == col)
		d =(*getrow(row)).dValue;
	else
	{
		cRowIterator conn = getrow(row); ++conn; // skip diag
		for(; !conn.isEnd(); ++conn)
		{
			if((*conn).iIndex == col)
			{
				d = (*conn).dValue;
				return;
			}		
		}	
		d = 0.0;
	}	
}
//! set whole matrix to a*I
//! used by Andreas to init matrix. avoid.
template<typename T>
bool SparseMatrix<T>::set(double a)
{
	for(int i=0; i<rows; i++)
		safeSetConnections(i, NULL);
	if(consmem)
	{
		delete[] consmem;
		consmem = NULL;
	}
	if(a == 0.0)
	{
		for(int i=0; i<rows; i++)
			iNrOfConnections[i] = 0;
		consmemsize = 0;
		iFragmentedMem = 0;
	}
	else
	{
		consmem = new connection[rows];
		consmemsize = rows;
		iFragmentedMem = 0;

		for(int i=0; i<rows; i++)
		{
			cons[i] = consmem+i;
			iNrOfConnections[i] = 1;
		}
	}
	return true;
}

#pragma mark private functions
////////////////////////////////////////////////////////////////////////////////


template<typename T>
void sortConnections(typename SparseMatrix<T>::connection *c, int nr, int row)
{
	if(nr<=1) return;
	// search diag
	if(c[0].iIndex != row)
		for(int i=1; i<nr; i++)
		{
			if(c[i].iIndex == row)
			{
				swap(c[i], c[0]);
				break;
			}			
		}
	if(nr-1 > 0)
		sort(c+1, c+1+nr-1);
}


#pragma mark row functions
////////////////////////////////////////////////////////////////////////////////


template<typename T>
void SparseMatrix<T>::setMatrixRow(int row, connection *c, int nr)
{
	sortConnections<entry_type>(c, nr, row);
	connection *n;
	if(c[0].iIndex != row)
	{
		n = new connection[nr+1];		
		n[0].iIndex = row; n[0].dValue = 0.0;
		for(int i=0; i<nr; i++)
			n[i+1] = c[i];
		nr++;
	}
	else
	{
		n = new connection[nr];
		for(int i=0; i<nr; i++)
			n[i] = c[i];
	}
	for(int i=0; i<nr; i++)
	{
		bandwidth = max(bandwidth, abs(n[i].iIndex - row));
		UG_ASSERT(n[i].iIndex >= 0 && n[i].iIndex < getRows(), *this << " cannot have connection " << n[i] << ".");
	}
	iFragmentedMem += nr;
	
	safeSetConnections(row, n);
	iTotalNrOfConnections += nr - iNrOfConnections[row];
	iNrOfConnections[row] = nr;
	iMaxNrOfConnections[row] = nr;
}

template<typename T>
void SparseMatrix<T>::addMatrixRow(int row, connection *c, int nr)
{	
	connection *old = cons[row];	
	if(old == NULL)
	{		
		setMatrixRow(row, c, nr);
		return;
	}
	
	if(nr == 1)
	{
		if(row == c->iIndex)
		{
			old[0].dValue += c->dValue;
			return;
		}

		int s;
		for(s=1; s<iNrOfConnections[row]; s++)
		{
			if(c->iIndex < old[s].iIndex) break;
		}

		if(c->iIndex == old[s].iIndex)
			old[s].dValue += c->dValue;
		else
		{
			int oldNrOfConnections = iNrOfConnections[row];
			connection *n = new connection[oldNrOfConnections+1];
			memcpy(n, old, sizeof(connection)*s);
			n[s] = *c;
			memcpy(n+s+1, old+s, sizeof(connection)*(oldNrOfConnections-s));
			
			safeSetConnections(row, n);
			iTotalNrOfConnections ++;
			iNrOfConnections[row]++;
			iMaxNrOfConnections[row] = iNrOfConnections[row];
		}
		return;
	}	
	
	//IFDEBUG( for(int i=0; i<nr; i++) UG_ASSERT(c[i].iIndex >= 0 && c[i].iIndex < getRows(), *A << " cannot have connection " << c[i] << "."); )

	int oldNrOfConnections = iNrOfConnections[row];
	// sort the connections
	sortConnections<entry_type>(c, nr, row);	
	
	int ic, skipped=0, iold=1;
	if(c[0].iIndex == row) 
	{
		old[0].dValue += c[0].dValue;
		ic = 1; 		
	}
	else
		ic = 0;
	
	while(ic < nr && iold < oldNrOfConnections)
	{
		if(c[ic].iIndex < old[iold].iIndex)
		{
			skipped++;
			ic++;
		}
		else if(c[ic].iIndex > old[iold].iIndex)
			iold++;
		else // "="
		{
			old[iold].dValue += c[ic].dValue;
			ic++;
			iold++;
		}
	}
	skipped += nr - ic;
	if(skipped == 0)  // everything already done
		return;
	

	// else realloc
	
	int iNewSize = oldNrOfConnections + skipped;
	connection *n = new connection[iNewSize];
	iFragmentedMem += iNewSize;
	n[0] = old[0];
	int i=1;
	iold=1;
	if(c[0].iIndex == row) 
		ic = 1; 		
	else
		ic = 0;

	while(ic < nr || iold < oldNrOfConnections)
	{
		if(iold >= oldNrOfConnections || (ic < nr && c[ic].iIndex < old[iold].iIndex))
		{
			UG_ASSERT(i < iNewSize && ic < nr, "");
			n[i++] = c[ic];
			ic++;
		}
		else if(ic >= nr || c[ic].iIndex > old[iold].iIndex)
		{
			UG_ASSERT(i < iNewSize && iold < oldNrOfConnections, "");
			n[i++] = old[iold];
			iold++;
		}
		else // "="
		{
			UG_ASSERT(i < iNewSize && iold < oldNrOfConnections, "");
			n[i++] = old[iold];
			ic++;
			iold++;
		}
	}

	safeSetConnections(row, n);
	iTotalNrOfConnections += iNewSize - oldNrOfConnections;
	iNrOfConnections[row] = iNewSize;
	iMaxNrOfConnections[row] = iNewSize;
}

template<typename T>
void SparseMatrix<T>::removezeros(int row)
{
	connection* con = cons[row];
	int nr = iNrOfConnections[row];
	// search from back for zero entries
	while(nr > 0 && con[nr-1].dValue == 0) nr--; // diagonaleintrag behalten
	for(int i=1; i < nr; i++)
		if(con[i].dValue == 0)
		{
			nr--;
			con[i] = con[nr];
			while(nr > 0 && con[nr-1].dValue == 0) nr--;
		}
	
	iNrOfConnections[row] = nr+1;
	iMaxNrOfConnections[row] = nr+1;
}


#pragma mark output functions
////////////////////////////////////////////////////////////////////////////////

#define CONNECTION_VIEWER_VERSION 1

// writeToFile
//--------------------------------------------------
//! writes to a file in somewhat SparseMatrix-market format (for connection viewer)
/*template<typename T>
void SparseMatrix<T>::writeToFile(const char *filename) const
{
	fstream file(filename, ios::out);
	file << CONNECTION_VIEWER_VERSION << endl;
	file << 2 << endl; // TODO: enter correct dimension here
	
	writePosToStream(file);	
	file << 1 << endl;
	for(int i=0; i < rows; i++)
	{
		for(cRowIterator conn = beginRow(i); !conn.isEnd(); ++conn)
			if((*conn).dValue != 0.0)
				file << GetOriginalIndex(tolevel, i) << " " << GetOriginalIndex(fromlevel, (*conn).iIndex) << " " << ((*conn).dValue) <<		endl;
	}
}*/

// writeToFile: in somewhat SparseMatrix-market format:
// length \n then for each connection: from to value
/*template<typename T>
 void SparseMatrix<T>::writeToFile2(const char *filename) const
 {
 fstream file(filename, ios::out);
 
 int level = min(fromlevel, tolevel);
 int m = max(rows, cols);
 cout << m << endl;
 for(int i=0; i < m; i++)
 {		
 int index = GetOriginalIndex(level, i)
 postype p = GetPosForIndex(index);
 file << p.x << " " << p.y << endl;
 }
 file << 1 << endl;
 for(int i=0; i < rows; i++)
 {
 for(cRowIterator conn = beginRow(i); !conn.isEnd(); ++conn)
 if((*conn).dValue != 0.0)
 file << GetOriginalIndex(tolevel, i) << " " << GetOriginalIndex(fromlevel, (*conn).iIndex) << " " << ((*conn).dValue) <<		endl;
 }
 }	*/


//!
//! print to console whole SparseMatrix
template<typename T>
void SparseMatrix<T>::print(const char * const text) const
{
	if(name) cout << endl << "================ " << name;
	if(text) cout << " == " << text;
	cout << " == fromlevel: " << fromlevel << " tolevel: " << tolevel << ", " << rows << "x" << cols << " =================" << endl;
	for(int i=0; i < rows; i++)
		getrow(i).print();				
}


//!
//! print the row row to the console
template<typename T>
void SparseMatrix<T>::printrow(int row) const
{
	cout << row << " [" << GetOriginalIndex(tolevel, row) << "]: ";
	for(int i=0; i < iNrOfConnections[row]; i++)
	{
		connection &c = cons[row][i];
		if(c.dValue == 0.0) continue;
		cout << " ";
		cout << "(" << c.iIndex << "[" << GetOriginalIndex(fromlevel, c.iIndex) << "]-> " << c.dValue << ")";			
	}
	//cout << " SUM: " << sum() << endl;
	cout << endl;
}

//! shortcut for gdb
template<typename T>
void SparseMatrix<T>::p() const
{
	print(NULL);
}

//! shortcut for gdb
template<typename T>
void SparseMatrix<T>::pr(int row) const
{
	printrow(row);
}

template<typename T>
void SparseMatrix<T>::printtype() const
{ 
	cout << *this; 
}






// safeSetConnections
//--------------------
//! "safe" way to set a connection, since when cons[row] is in the big consecutive consmem-array,
//! you mustn'd delete it.
template<typename T>
void SparseMatrix<T>::safeSetConnections(int row, connection *mem) const
{	
	if(cons[row] != NULL && cons[row] < consmem || cons[row] > consmem + consmemsize)
		delete[] cons[row];
	cons[row] = mem;
}


#pragma mark Template Expressions
////////////////////////////////////////////////////////////////////////////////

// operator *
//--------------------
//!
//! * operator for template expression x = A y. (x, y Vector, A SparseMatrix).
/*
template<typename entry_type, typename Vector_type>
Expression<SparseMatrix<entry_type>, Multiply_Operator<entry_type, typename Vector_type::entry_type>, Vector_type> 
	operator * (const SparseMatrix<entry_type> &l, const XD< Vector_type > &r)
{ 
	return Expression<SparseMatrix<entry_type>, Multiply_Operator<entry_type, typename Vector_type::entry_type>, Vector_type> (l, r.cast()); 
}


// todo: prevent x = A * x; mit feld forbiddenDestination
//!
//! class for template Expression x = A y (x, y Vector, A SparseMatrix).
template<typename entry_type, typename Vector_type> 
class Expression<SparseMatrix<entry_type>, Multiply_Operator<entry_type, typename Vector_type::entry_type >, Vector_type > 
: public XD< Expression<SparseMatrix<entry_type>, Multiply_Operator<entry_type, typename Vector_type::entry_type >, Vector_type > >
{ 
public:	
	typedef typename Multiply_Operator<entry_type, typename Vector_type::entry_type>::ReturnType ReturnType;
	const SparseMatrix<entry_type>& l;
	const Vector_type & r; 
	inline Expression(const SparseMatrix<entry_type> &l_, const Vector_type &r_) : l(l_), r(r_) {} 
	
	inline ReturnType operator [] (int i) const
	{
		return l[i] * r;
	} 
	
	inline void copyTo(ReturnType &d, int i) const
	{
		l[i].assign_mult(d, r);
	}
	
	inline void addTo(ReturnType &d, int i) const
	{
		l[i].add_mult(d, r);
	}
	
	inline void substractFrom(ReturnType &d, int i) const
	{
		l[i].sub_mult(d, r);
	}	
	
	inline int size() const	{	return l.size();	}
	
	// print routines
	friend ostream &operator<<(ostream &output, 
			const Expression<SparseMatrix<entry_type>, Multiply_Operator<entry_type, typename Vector_type::entry_type >, Vector_type >  &ex)
	{
		output << "(" << ex.l	<< "*" << ex.r << ")"; 
		return output;
	}
	inline void printtype() const	{	cout << *this; }
}; 

/*
// todo: prevent x = A * x; mit feld forbiddenDestination
// x = r / A.Diag();
/*
template<> 
struct Expression<Vector, Divide_Operator, SparseMatrix<T>::diagcomponent> 
{ 
	const Vector& l;
	const SparseMatrix<T>::diagcomponent& r; 
	inline Expression(const Vector &l_, const SparseMatrix<T>::diagcomponent &r_) : l(l_), r(r_) 
	{ UG_ASSERT(l.size() == r.size(), l << " has different length as " <<  r); }
	
	inline double operator [] (int i) const
	{
		return l[i] / r[i];
	} 
	
	inline int size() const	{	return l.size();	}
	
	// print routines
	friend ostream &operator<<(ostream &output, const Expression<Vector, Divide_Operator, SparseMatrix<T>::diagcomponent>  &ex)
	{
		output << "(" << ex.l << Divide_Operator::cTyp() << ex.r << ")"; 
		return output;
	}
	inline void printtype() const	{	cout << *this; }
}; 

Expression<Vector, Divide_Operator, SparseMatrix<T>::diagcomponent> operator/(const Vector &l, const SparseMatrix<T>::diagcomponent &r);
*/	



//!
//!
//! defrag
template<typename T>
void SparseMatrix<T>::defrag()
{
	//UG_ASSERT(0, "this function is broken");
	iTotalNrOfConnections=0;
	for(int i=0; i<rows; i++)
		iTotalNrOfConnections+=iNrOfConnections[i];
	
	connection *consmemNew = new connection[iTotalNrOfConnections+3];
	connection *p= consmemNew;
	for(int i=0; i<rows; i++)
	{
		for(int k=0; k<iNrOfConnections[i]; k++)
			swap(p[k], cons[i][k]);
		safeSetConnections(i, p);
		p += iNrOfConnections[i];
		iMaxNrOfConnections[i] = iNrOfConnections[i];
	}
	delete[] consmem;
	consmem = consmemNew;
	iFragmentedMem = 0;
	consmemsize = iTotalNrOfConnections;
}


//!
//! @param node
//! @param depth 
template<typename T>
void SparseMatrix<T>::getNeighborhood(int node, int depth, vector<int> &indices, int *posInConnections) const
{
	// do this with a map???
	indices.clear();
	
	vector<cRowIterator> iterators;
	iterators.reserve(depth);
	
	iterators.push_back( beginRow(node) );
	
	while(iterators.size() != 0)
	{
		if(iterators.back().isEnd())
			iterators.pop_back();			
		else
		{
			int index = (*iterators.back()).iIndex;
			++iterators.back();
			if(iterators.size() < depth)
				iterators.push_back( beginRow(index) );
			else
			{
				int pos;
				if(posInConnections == NULL)
				{
					for(pos=0; pos<indices.size(); pos++)
						if(indices[pos] == index)
							break;
					if(pos == indices.size())
						indices.push_back(index);
				}
				else 
				{
					pos = posInConnections[index];
					if(pos == -1)
					{
						pos = posInConnections[index] = indices.size();
						indices.push_back(index);
					}
				}
				// else (count etc.)
			}				
									
			
		}
	}
	
	if(posInConnections)
	{
		for(int i=0; i<indices.size(); i++)
			posInConnections[indices[i]] = -1;
	}	

	sort(indices.begin(), indices.end());

}




template<typename T>
bool SparseMatrix<T>::isCloseToBoundary(int node, int distance) const
{
	if(distance == 0) return isUnconnected(node);
	bool bFound = false;
	for(cRowIterator itA = beginRow(node); !itA.isEnd() && !bFound; ++itA)
		bFound = isCloseToBoundary((*itA).iIndex, distance-1);
		
	return bFound;	
}

} // namespace ug



