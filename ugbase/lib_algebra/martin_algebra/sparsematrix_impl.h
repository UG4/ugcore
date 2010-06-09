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
#include "algebra_misc.h"

template<typename T> T abs(const T &a, const T &b)
{
	if(a > b) return a-b; else return b-a;
}

////////////////////////////////////////////////////////////////////////////////
static size_t GetOriginalIndex(size_t level, size_t i) { return i; }

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
	for(size_t i=0; i < rows; i++)
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
bool SparseMatrix<T>::create(size_t _rows, size_t _cols)
{
	UG_ASSERT(rows == 0 && cols == 0, *this << " not empty.");
	
	rows = _rows;
	cols = _cols;
	
	cons = new connection*[rows+1];
	memset(cons, 0, sizeof(connection*)*(rows+1));

	iNrOfConnections = new size_t[rows];
	memset(iNrOfConnections, 0, sizeof(size_t)*rows);
	
	iMaxNrOfConnections = new size_t[rows];
	memset(iMaxNrOfConnections, 0, sizeof(size_t)*rows);
	
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
	for(size_t j=0; j < B.rows; j++)
		for(cRowIterator conn = B.beginRow(j); !conn.isEnd(); ++conn)
			if((*conn).dValue == 0) continue;
			else iNrOfConnections[(*conn).iIndex]++;
	
	size_t newTotal = 0;
	for(size_t i=0; i < rows; i++)
		newTotal += iNrOfConnections[i];
	
	
	
	size_t *nr = new size_t[rows];
	// init SparseMatrix data structure
	consmem = new connection[newTotal];
	consmemsize = newTotal;
	
	connection *p = consmem;
	for(size_t i=0; i < rows; i++)
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
	for(size_t i=0; i < B.rows; i++)
		for(cRowIterator conn = B.beginRow(i); !conn.isEnd(); ++conn)
		{
			if((*conn).dValue == 0) continue;
			size_t ndx = (*conn).iIndex;
			UG_ASSERT(ndx >= 0 && ndx < rows, "connection " << (*conn) << " of " << B << ", row " << i << " out of range 0.." << rows);
			
			size_t k= nr[ndx];
			
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






/*void SparseMatrix<T>::recreateWithMaxNrOfConnections(size_t newMax) const
 {
 // create new cons Memory
 connection *consmemNew = new connection[newMax];
 
 // adjust pointers
 size_t diff = consmemNew - consmem;
 for(size_t i=0; i<rows; i++)
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
	for(size_t i=0; i<rows; i++)
	{
		if(isUnconnected(i)) continue;		
		for(rowIterator conn = beginRow(i); !conn.isEnd(); ++conn)
		{
			size_t conindex = (*conn).iIndex;
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
void SparseMatrix<T>::setDirichletRow(size_t row)
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
void SparseMatrix<T>::setDirichletRows(size_t *pRows, size_t nrows)
{
	for(size_t i=0; i<nrows; i++)
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
	size_t n = min(rows, cols);
	for(size_t i=0; i < n; i++)
	{
		cRowIterator conn = beginRow(i);
		res[i] = 0.0;
		for(; !conn.isEnd(); ++conn)
			res[i] += (*conn).dValue * x[(*conn).iIndex];
	}
	
	for(size_t i=n; n < cols; i++)
	{
		cRowIterator conn = beginRow(i);
		++conn; // skip diag, since x.size <= n
		res[i] = 0.0;
		for(; !conn.isEnd(); ++conn)
			res[i] += (*conn).dValue * x[(*conn).iIndex];
	}*/

	for(size_t i=0; i < rows; i++)
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
	size_t n = min(rows, x.size());
	for(size_t i=0; i<n; i++)
	{
		cRowIterator conn = beginRow(i);
		for(; !conn.isEnd(); ++conn)
			res[(*conn).iIndex] += (*conn).dValue * x[i];
	}
	*/

	/*size_t n = min(cols, rows);
	for(size_t i=0; i<n; i++)
	{
		cRowIterator conn = beginRow(i);
		for(; !conn.isEnd(); ++conn)
			res[(*conn).iIndex] += (*conn).dValue * x[i];
	}

	for(size_t i=n; i<rows; i++)
	{
		cRowIterator conn = beginRow(i);
		++conn;
		for(; !conn.isEnd(); ++conn)
			res[(*conn).iIndex] += (*conn).dValue * x[i];
	}*/
	for(size_t i=0; i<rows; i++)
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
	
	size_t n = min(rows, x.size());
	for(size_t i=0; i < n; i++)
	{
		cRowIterator conn = beginRow(i);
		for(; !conn.isEnd(); ++conn)
			res[i] -= (*conn).dValue * x[(*conn).iIndex];
	}
	
	for(size_t i=n; n < rows; i++)
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
const matrixrow<T> SparseMatrix<T>::getrow(size_t i) const
{	
	return matrixrow<entry_type> (*this, i);
}

//!
//! @return a matrixrow object of row i
template<typename T>
const matrixrow<T> SparseMatrix<T>::operator [] (size_t i) const
{		
	return matrixrow<entry_type> (*this, i);
}


// getDiag
//-------------
//! get Diagonal A_[i,i] of matrix
//! @param i
//! @return A_{i,i}
template<typename T>
inline const T &SparseMatrix<T>::getDiag(size_t i) const
{
	UG_ASSERT(cons[i][0].iIndex == i, *this << " first entry has to be diagonal");
	// evtl anders, da nicht jede Matrix diageintrag
	return cons[i][0].dValue;
}

template<typename T>
inline T &SparseMatrix<T>::getDiag(size_t i)
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
inline bool SparseMatrix<T>::isUnconnected(size_t i) const
{
	UG_ASSERT(i < rows && i >= 0, *this << ": " << i << " out of bounds.");
	return iNrOfConnections[i] <= 1; 
}




// add Submatrix
//--------------------
//! function to add submatrices ( submatrix )
template<typename T>
template<typename M>
void SparseMatrix<T>::add(const M &mat, size_t *rows, size_t *cols)
{
	connection *c = new connection[mat.getCols()];
	size_t nc;
	for(size_t i=0; i < mat.getRows(); i++)
	{
		nc = 0;
		for(size_t j=0; j < mat.getCols(); j++)
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
void SparseMatrix<T>::set(const M &mat, size_t *rows, size_t *cols)
{
	connection *c = new connection[mat.getCols()];
	size_t nc;
	for(size_t i=0; i < mat.getRows(); i++)
	{
		nc = 0;
		for(size_t j=0; j < mat.getCols(); j++)
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
void SparseMatrix<T>::get(M &mat, size_t *rows, size_t *cols) const
{
	vector<sortStruct<size_t> > sortedCols(mat.getCols());
	
	for(size_t i=0; i<mat.getCols(); i++)
	{
		sortedCols[i].index = i;
		sortedCols[i].sortValue = cols[i];
	}	
	sort(sortedCols.begin(), sortedCols.end());
	
	for(size_t i=0; i < mat.getRows(); i++)
	{
		size_t iindex_global = rows[i];
		size_t iindex_local = i;
	
		cRowIterator conn = beginRow(iindex_global);
		// diagonal
		mat(iindex_local, iindex_local) = (*conn).dValue;
		++conn;

		size_t j = 0;
		while(j < mat.getCols() && !conn.isEnd())
		{
			size_t jindex_global = sortedCols[j].sortValue;
			size_t jindex_local = sortedCols[j].index;

			if(iindex_global == jindex_global) { j++; continue; }
			size_t cindex_global = (*conn).iIndex;
			
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
void SparseMatrix<T>::add(const M &mat, vector<size_t> &rows, vector<size_t> &cols)
{
	ASSERT1(mat.getCols() == cols.size() && mat.getRows() == rows.size());
	add(mat, &rows[0], &cols[0]);
}

template<typename T>
template<typename M>
void SparseMatrix<T>::set(const M &mat, vector<size_t> &rows, vector<size_t> &cols)
{
	ASSERT1(mat.getCols() == cols.size() && mat.getRows() == rows.size());
	set(mat, &rows[0], &cols[0]);
}

template<typename T>
template<typename M>
void SparseMatrix<T>::get(M &mat, vector<size_t> &rows, vector<size_t> &cols) const
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
	for(size_t i=0; i < I.size(); i++)
	{
		nc = 0;
		for(size_t j=0; j < J.size(); j++)
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
void SparseMatrix<T>::add(const entry_type &d, size_t row, size_t col)
{
	connection c;
	c.dValue = d;
	c.iIndex = col;
	addMatrixRow(row, &c, 1);
}

template<typename T>
void SparseMatrix<T>::set(const entry_type &d, size_t row, size_t col)
{
	connection c;
	c.dValue = d;
	c.iIndex = col;
	setMatrixRow(row, &c, 1);
}

template<typename T>
void SparseMatrix<T>::get(entry_type &d, size_t row, size_t col) const
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
	for(size_t i=0; i<rows; i++)
		safeSetConnections(i, NULL);
	if(consmem)
	{
		delete[] consmem;
		consmem = NULL;
	}
	if(a == 0.0)
	{
		for(size_t i=0; i<rows; i++)
			iNrOfConnections[i] = 0;
		consmemsize = 0;
		iFragmentedMem = 0;
	}
	else
	{
		consmem = new connection[rows];
		consmemsize = rows;
		iFragmentedMem = 0;

		for(size_t i=0; i<rows; i++)
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
void sortConnections(typename SparseMatrix<T>::connection *c, size_t nr, size_t row)
{
	if(nr<=1) return;
	// search diag
	if(c[0].iIndex != row)
		for(size_t i=1; i<nr; i++)
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
void SparseMatrix<T>::setMatrixRow(size_t row, connection *c, size_t nr)
{
	sortConnections<entry_type>(c, nr, row);
	connection *n;
	if(c[0].iIndex != row)
	{
		n = new connection[nr+1];		
		n[0].iIndex = row; n[0].dValue = 0.0;
		for(size_t i=0; i<nr; i++)
			n[i+1] = c[i];
		nr++;
	}
	else
	{
		n = new connection[nr];
		for(size_t i=0; i<nr; i++)
			n[i] = c[i];
	}
	for(size_t i=0; i<nr; i++)
	{
		bandwidth = max(bandwidth, abs(n[i].iIndex, row));
		UG_ASSERT(n[i].iIndex >= 0 && n[i].iIndex < getRows(), *this << " cannot have connection " << n[i] << ".");
	}
	iFragmentedMem += nr;
	
	safeSetConnections(row, n);
	iTotalNrOfConnections += nr - iNrOfConnections[row];
	iNrOfConnections[row] = nr;
	iMaxNrOfConnections[row] = nr;
}

template<typename T>
void SparseMatrix<T>::addMatrixRow(size_t row, connection *c, size_t nr)
{	
	if(nr <= 0) return;

	//cout << "row: " << row << " nr: " << nr << endl;
	connection *old = cons[row];	
	if(old == NULL)
	{		
		setMatrixRow(row, c, nr);
		return;
	}
	
	UG_ASSERT(iNrOfConnections[row] != 0, "cons[row] != NULL but iNrOfConnections[row] == 0 ???");

	// the matrix row is not empty and we are adding more than one connection

	size_t oldNrOfConnections = iNrOfConnections[row];

	// sort the connections: diagonal, then index1 < index2 < ...
	sortConnections<entry_type>(c, nr, row);	
	
	// diagonal
	size_t ic, skipped=0, iold=1;
	if(c[0].iIndex == row) 
	{
		old[0].dValue += c[0].dValue;
		ic = 1; 		
	}
	else
		ic = 0;
	
	// off-diagonals:
	// - add connections which are there
	// - count not existing connections (=skipped)
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
	// add the rest (when iold == oldNrOfConnections, but still ic < nr)
	skipped += nr - ic;

	if(skipped == 0)  // everything already done
		return;
	
	// else realloc
	
	size_t iNewSize = oldNrOfConnections + skipped;
	connection *n = new connection[iNewSize];
	iFragmentedMem += iNewSize;

	size_t i=0; iold=0;
	// deal with diagonal
	n[i++] = old[iold++];
	if(c[0].iIndex == row) 
		ic = 1; 		
	else
		ic = 0;

	// merge the two arrays
	while(ic < nr && iold < oldNrOfConnections)
	{
		if(c[ic].iIndex < old[iold].iIndex)
		{
			UG_ASSERT(c[ic].iIndex >= 0 && c[ic].iIndex < getCols(), "");
			n[i++] = c[ic++];
		}
		else if(c[ic].iIndex > old[iold].iIndex)
			n[i++] = old[iold++];
		else // "="
		{
			// "add" already done
			n[i++] = old[iold];
			ic++;
			iold++;
		}
	}
	// deal with the rest
	while(ic<nr) n[i++] = c[ic++];
	while(iold<oldNrOfConnections) n[i++] = old[iold++];

	UG_ASSERT(i == iNewSize, "row: " << row << " i: " << i << " iNewSize: " << iNewSize);

	// set connections
	safeSetConnections(row, n);
	iTotalNrOfConnections += iNewSize - oldNrOfConnections;
	iNrOfConnections[row] = iNewSize;
	iMaxNrOfConnections[row] = iNewSize;
}

template<typename T>
void SparseMatrix<T>::removezeros(size_t row)
{
	connection* con = cons[row];
	size_t nr = iNrOfConnections[row];
	// search from back for zero entries
	while(nr > 0 && con[nr-1].dValue == 0) nr--; // diagonaleintrag behalten
	for(size_t i=1; i < nr; i++)
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


//!
//! print to console whole SparseMatrix
template<typename T>
void SparseMatrix<T>::print(const char * const text) const
{
	if(name) cout << endl << "================ " << name;
	if(text) cout << " == " << text;
	cout << " == fromlevel: " << fromlevel << " tolevel: " << tolevel << ", " << rows << "x" << cols << " =================" << endl;
	for(size_t i=0; i < rows; i++)
		getrow(i).print();				
}


//!
//! print the row row to the console
template<typename T>
void SparseMatrix<T>::printrow(size_t row) const
{
	cout << row << " [" << GetOriginalIndex(tolevel, row) << "]: ";
	for(size_t i=0; i < iNrOfConnections[row]; i++)
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
void SparseMatrix<T>::pr(size_t row) const
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
void SparseMatrix<T>::safeSetConnections(size_t row, connection *mem) const
{	
	if(cons[row] != NULL && cons[row] < consmem || cons[row] > consmem + consmemsize)
		delete[] cons[row];
	cons[row] = mem;
}


#pragma mark Template Expressions
////////////////////////////////////////////////////////////////////////////////



//!
//!
//! defrag
template<typename T>
void SparseMatrix<T>::defrag()
{
	iTotalNrOfConnections=0;
	for(size_t i=0; i<rows; i++)
		iTotalNrOfConnections+=iNrOfConnections[i];
	
	connection *consmemNew = new connection[iTotalNrOfConnections+3];
	connection *p= consmemNew;
	for(size_t i=0; i<rows; i++)
	{
		for(size_t k=0; k<iNrOfConnections[i]; k++)
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




} // namespace ug



