/*
 *  SparseMatrix.hpp
 *  flexamg
 *
 *  Created by Martin Rupp on 04.11.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 */
#pragma once
// creation etc

#include <fstream>

#ifndef FLEXAMG
#include "algebra_misc.h"
#endif

template<typename T> T abs(const T &a, const T &b)
{
	if(a > b) return a-b; else return b-a;
}


namespace ug{
//!
//! constructor for empty SparseMatrix
template<typename T>
SparseMatrix<T>::SparseMatrix()
{
	name = "?";
	cols = rows = iTotalNrOfConnections = 0;
	iFragmentedMem = 0;
	pRowStart = pRowEnd = NULL;
	consmem = NULL; consmemsize = 0;
	iTotalNrOfConnections = 0;
	bandwidth = 0;
	
	estimatedRowSize = 0;
	iMaxNrOfConnections = NULL;
	
#ifdef FLEXAMG
	tolevel = fromlevel = 0;
#endif

	FORCE_CREATION { print(); p(); pr(0); }
}

template<typename T>
bool SparseMatrix<T>::isFinalized() const
{
	return pRowEnd == (pRowStart+1);
}

template<typename T> bool
SparseMatrix<T>::destroy()
{
	if(isFinalized())
	{
		if(pRowStart) delete [] pRowStart;
		pRowStart = NULL;
		pRowEnd = NULL;
	}
	else
	{
		if(pRowStart)
		{
			for(size_t i=0; i < rows; i++)
				safeSetConnections(i, NULL);
			delete [] pRowStart;
			pRowStart = NULL;
		}
		if(pRowEnd) delete [] pRowStart; pRowStart = NULL;
	}

	if(consmem)	delete [] consmem; consmem = NULL;
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
	
	pRowStart = new connection*[rows+1];
	memset(pRowStart, 0, sizeof(connection*)*(rows+1));

	pRowEnd = new connection*[rows+1];
	memset(pRowEnd, 0, sizeof(connection*)*(rows+1));
	

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
	UG_ASSERT(rows == 0 && cols == 0, *this << " not empty.");
	
	rows = B.cols;
	cols = B.rows;
#ifdef FLEXAMG
	tolevel = B.fromlevel;
	fromlevel = B.tolevel;
#endif

	pRowStart = new connection*[rows+1];
	pRowEnd= pRowStart+1;
	iMaxNrOfConnections = new size_t[rows];
	iTotalNrOfConnections = 0;
	bandwidth = 0;

	//-----------------------------------

	size_t *nr = new size_t[rows];
	memset(nr, 0, sizeof(size_t)*rows);

	// get length of each row
	for(size_t j=0; j < B.rows; j++)
		for(cRowIterator conn = B.beginRow(j); !conn.isEnd(); ++conn)
			if((*conn).dValue == 0) continue;
			else nr[(*conn).iIndex]++;
	

	size_t newTotal = 0;
	for(size_t i=0; i < rows; i++)
		newTotal += nr[i];
	
	
	// allocate one big connection memory
	consmem = new connection[newTotal];
	consmemsize = newTotal;
	
	// init pRowStart array
	// (this also inits pRowEnd)
	connection *p = consmem;
	for(size_t i=0; i < rows; i++)
	{
		pRowStart[i] = p;
		p += nr[i];
		iMaxNrOfConnections[i] = nr[i];
	}
	pRowStart[rows] = p;


	// write values
	memset(nr, 0, sizeof(size_t)*rows);
	iTotalNrOfConnections = newTotal;
	iFragmentedMem = 0;
	
	bandwidth = 0;

	for(size_t i=0; i < B.rows; i++)
		for(cRowIterator conn = B.beginRow(i); !conn.isEnd(); ++conn)
		{
			if((*conn).dValue == 0) continue;
			size_t ndx = (*conn).iIndex;
			UG_ASSERT(ndx >= 0 && ndx < rows, "connection " << (*conn) << " of " << B << ", row " << i << " out of range 0.." << rows);
			
			size_t k = nr[ndx];
			
			UG_ASSERT(k>=0 && k<getNrOfConnections(ndx), "k = " << k << ". precalculated nr of Connections " << getNrOfConnections(ndx) << " wrong?");
			
			pRowStart[ndx][k].dValue = (*conn).dValue;
			pRowStart[ndx][k].iIndex = i;
			if(bandwidth < i-ndx) bandwidth = i-ndx;
			nr[ndx]++;
		}
	
	// note that all connections are sorted.

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



//general functions
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
				SubMult(b[i], (*conn).dValue, b[conindex]);
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
	if(getNrOfConnections(row) == 0 || !isInConsMem(row))
	{
		if(isInConsMem(row))
			iFragmentedMem += 1;
		else
			iFragmentedMem += 1-getNrOfConnections(row);

		connection *c = new connection[1];
		safeSetConnections(row, c);
		iMaxNrOfConnections[row] = 1;
	}

	pRowStart[row][0].dValue = 1.0;
	pRowStart[row][0].iIndex= row;
	pRowEnd[row] = pRowStart[row]+1;

}

//!
//! sets the # nrows in pRows to i,i = 1.0.
template<typename T>
void SparseMatrix<T>::setDirichletRows(size_t *pRows, size_t nrows)
{
	for(size_t i=0; i<nrows; i++)
		setDirichletRow(pRows[i]);
}

#ifndef FLEXAMG
template<typename T>
bool SparseMatrix<T>::set_dirichlet_rows(const local_index_type &ind)
{
	for(std::size_t i=0; i<ind.size(); i++)
		setDirichletRow(ind[i][0]);
	return true;
}
#endif


// res = A*x
template<typename T>
template<typename Vector_type>
bool SparseMatrix<T>::apply(Vector_type &res, const Vector_type &x) const
{
	UG_ASSERT(cols == x.size(), "x: " << x << " has wrong length (should be " << cols << "). A: " << *this);
	UG_ASSERT(rows == res.size(), "res: " << x << " has wrong length (should be " << rows << "). A: " << *this);
	
	for(size_t i=0; i < rows; i++)
	{
		res[i] = 0.0;
		for(cRowIterator conn = beginRow(i); !conn.isEnd(); ++conn)
			res[i] += (*conn).dValue * x[(*conn).iIndex];
	}

	return true;
}



// res = A.T() * x
template<typename T>
template<typename Vector_type>
bool SparseMatrix<T>::apply_transposed(Vector_type &res, const Vector_type &x) const
{
	UG_ASSERT(rows == x.size(), "x: " << x << " has wrong length (should be " << rows << "). A: " << *this);
	UG_ASSERT(cols == res.size(), "res: " << x << " has wrong length (should be " << cols << "). A: " << *this);
	
	res = 0.0;

	for(size_t i=0; i<rows; i++)
	{
		for(cRowIterator conn = beginRow(i); !conn.isEnd(); ++conn)
		{
			if((*conn).dValue != 0.0)
				res[(*conn).iIndex] += (*conn).dValue * x[i];
		}
	}
	return true;
}


// res = res - A*x
template<typename T>
template<typename Vector_type>
bool SparseMatrix<T>::matmul_minus(Vector_type &res, const Vector_type &x) const
{
	UG_ASSERT(cols == x.size(), "x: " << x << " has wrong length (should be " << cols << "). A: " << *this);
	UG_ASSERT(rows == res.size(), "res: " << x << " has wrong length (should be " << rows << "). A: " << *this);
	
	for(size_t i=0; i < rows; i++)
	{
		cRowIterator conn = beginRow(i);
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
	return operator() (i, i);
}

template<typename T>
inline T &SparseMatrix<T>::getDiag(size_t i)
{
	return operator() (i, i);
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
	if(pRowStart[i] == NULL)
		return true;
	int nr=getNrOfConnections(i);
	if(nr == 0 || (pRowStart[i][0].iIndex == i && nr == 1))
		return true;

	return false;
}

template<typename T>
inline size_t SparseMatrix<T>::getNrOfConnections(size_t row) const
{
	return pRowEnd[row] - pRowStart[row];
}


//! set whole matrix to a*I
//! used by Andreas to init matrix. avoid.
template<typename T>
bool SparseMatrix<T>::set(double a)
{
	if(a == 0.0)
	{
		UG_ASSERT(pRowStart != NULL, "");
		for(size_t row=0; row<rows; row++)
			for(rowIterator it = beginRow(row); !it.isEnd(); ++it)
				(*it).dValue = 0.0;
	}
	else
	{
		UG_ASSERT(0, "not implemented fully");
		for(size_t row=0; row<rows; row++)
			for(rowIterator it = beginRow(row); !it.isEnd(); ++it)
			{
				if((*it).iIndex == row)
					(*it).dValue = a;
				else
					(*it).dValue = 0.0;
			}
	}
	return true;
}


////////////////////////////////////////////////////////////////////////////////


template<typename T>
void sortConnections(typename SparseMatrix<T>::connection *c, size_t nr, size_t row)
{

}


////////////////////////////////////////////////////////////////////////////////


template<typename T>
void SparseMatrix<T>::setMatrixRow(size_t row, connection *c, size_t nr)
{
	definalize();

	if(isInConsMem(row))
		iFragmentedMem += nr;
	else
		iFragmentedMem += nr - getNrOfConnections(row);

	connection *n = new connection[nr];
	for(size_t i=0; i<nr; i++)
		n[i] = c[i];

	sort(n, n+nr);

	for(size_t i=0; i<nr; i++)
		UG_ASSERT(n[i].iIndex >= 0 && n[i].iIndex < num_rows(), *this << " cannot have connection " << n[i] << ".");


	bandwidth = max(bandwidth, abs(n[0].iIndex, row));
	bandwidth = max(bandwidth, abs(n[nr-1].iIndex, row));
	

	safeSetConnections(row, n);
	iTotalNrOfConnections += nr - getNrOfConnections(row);

	pRowEnd[row] = pRowStart[row]+nr;
	iMaxNrOfConnections[row] = nr;
}

template<typename T>
void SparseMatrix<T>::addMatrixRow(size_t row, connection *c, size_t nr)
{	
	if(nr <= 0) return;

	//cout << "row: " << row << " nr: " << nr << endl;
	connection *old = pRowStart[row];

	if(old == NULL)
	{		
		setMatrixRow(row, c, nr);
		return;
	}
	
	UG_ASSERT(getNrOfConnections(row) != 0, "cons[row] != NULL but getNrOfConnections(row) == 0 ???");

	// the matrix row is not empty and we are adding more than one connection

	size_t oldNrOfConnections = getNrOfConnections(row);

	// sort the connections: diagonal, then index1 < index2 < ...
	sort(c, c+nr);
	
	size_t ic=0, iold=0, skipped=0;
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
	definalize();
	
	size_t iNewSize = oldNrOfConnections + skipped;
	connection *n = new connection[iNewSize];

	if(isInConsMem(row))
		iFragmentedMem += iNewSize;
	else
		iFragmentedMem += iNewSize - getNrOfConnections(row);

	size_t i=0; iold=0; ic = 0;

	// merge the two arrays
	while(ic < nr && iold < oldNrOfConnections)
	{
		if(c[ic].iIndex < old[iold].iIndex)
		{
			UG_ASSERT(c[ic].iIndex >= 0 && c[ic].iIndex < num_cols(), "");
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
	pRowEnd[row] = pRowStart[row]+iNewSize;
	iMaxNrOfConnections[row] = iNewSize;
}



////////////////////////////////////////////////////////////////////////////////


//!
//! print to console whole SparseMatrix
template<typename T>
void SparseMatrix<T>::print(const char * const text) const
{
	cout << "================= SparseMatrix " << rows << "x" << cols << " =================" << endl;
	for(size_t i=0; i < rows; i++)
		getrow(i).print();				
}


//!
//! print the row row to the console
template<typename T>
void SparseMatrix<T>::printrow(size_t row) const
{
#ifdef FLEXAMG
	cout << "row " << row << " [" << GetOriginalIndex(tolevel, row) << "] : ";
#else
	cout << "row " << row << ": ";
#endif
	for(cRowIterator it=beginRow(row); !it.isEnd(); ++it)
	{
		if((*it).dValue == 0.0) continue;
		cout << " ";
#ifdef FLEXAMG
		cout << "(" << (*it).iIndex << " [" << GetOriginalIndex(fromlevel, (*it).iIndex) << "] -> " << (*it).dValue << ")";
#else
		cout << "(" << (*it).iIndex << " -> " << (*it).dValue << ")";
#endif
	}
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


#define CONNECTION_VIEWER_VERSION 1

#ifdef FLEXAMG
// writeToFile
//--------------------------------------------------
//! writes to a file in somewhat SparseMatrix-market format (for connection viewer)
template<typename T>
void SparseMatrix<T>::writeToFile(const char *filename) const
{
	fstream file(filename, ios::out);
	file << CONNECTION_VIEWER_VERSION << endl;
	file << flexamg_dimensions << endl;

	writePosToStream(file);
	file << 1 << endl;
	for(int i=0; i < rows; i++)
	{
		for(cRowIterator conn = beginRow(i); !conn.isEnd(); ++conn)
			if((*conn).dValue != 0.0)
				file << GetOriginalIndex(tolevel, i) << " " << GetOriginalIndex(fromlevel, (*conn).iIndex) << " " << ((*conn).dValue) <<		endl;
	}
}

#endif


////////////////////////////////////////////////////////////////////////////////


template<typename T>
bool SparseMatrix<T>::isInConsMem(size_t row) const
{
	return pRowStart[row] >= consmem && pRowEnd[row] < consmem + consmemsize;
}


// safeSetConnections
//--------------------
//! "safe" way to set a connection, since when cons[row] is in the big consecutive consmem-array,
//! you mustn'd delete it.
template<typename T>
void SparseMatrix<T>::safeSetConnections(size_t row, connection *mem) const
{	
	if(pRowStart[row] != NULL && !isInConsMem(row))
		delete[] pRowStart[row];
	pRowStart[row] = mem;
}


template<typename T>
void SparseMatrix<T>::definalize()
{
	if(!isFinalized())
		return; // already not finalized.
	connection **p = new connection*[rows];
	memcpy(p, pRowEnd, sizeof(connection*)*rows);
	pRowEnd = p;
}

//!
//!
//! defrag
template<typename T>
void SparseMatrix<T>::defrag()
{
	iTotalNrOfConnections=0;
	for(size_t i=0; i<rows; i++)
		iTotalNrOfConnections+= getNrOfConnections(i);
	
	connection *consmemNew = new connection[iTotalNrOfConnections];
	connection *p = consmemNew;
	for(size_t i=0; i<rows; i++)
	{
		int nr=getNrOfConnections(i);
		for(size_t k=0; k < nr; k++)
			swap(p[k], pRowStart[i][k]);
		safeSetConnections(i, p);
		p += nr;
		iMaxNrOfConnections[i] = nr;
	}
	pRowStart[rows] = p;

	delete[] pRowEnd;
	pRowEnd = pRowStart+1;

	delete[] consmem;
	consmem = consmemNew;

	iFragmentedMem = 0;
	consmemsize = iTotalNrOfConnections;
}



/*template<typename T>
 void SparseMatrix<T>::removezeros(size_t row)
 {
 connection* con = cons[row];
 size_t nr = iNrOfConnections[row];
 // search from back for zero entries
 for(size_t i=1; i < nr; i++)
 if(con[i].iIndex != row && con[i].dValue == 0)
 {
 nr--;
 con[i] = con[nr];
 while(nr > 0 && con[nr-1].dValue == 0) nr--;
 }

 iNrOfConnections[row] = nr+1;
 iMaxNrOfConnections[row] = nr+1;
 }*/

template<typename T>
size_t SparseMatrix<T>::getConnection(size_t r, size_t c) const
{
	connection *con = pRowStart[r];
	int nr = getNrOfConnections(r);
	size_t left=0, right = nr-1;
	while(left <= right)
	{
		size_t mid = (left+right)/2;
		if (con[mid].iIndex < c)
			left = mid+1;
		else if (con[mid].iIndex > c)
			right = mid-1;
		else
			return mid;
	}
	return -1;
}

template<typename T>
const typename SparseMatrix<T>::entry_type &SparseMatrix<T>::operator() (int r, int c) const
{
	int i = getConnection(r, c);
	if(i==-1)
	{
		static T t=0;
		return t;
	}
	else
		return pRowStart[r][i].dValue;
}

template<typename T>
typename SparseMatrix<T>::entry_type &SparseMatrix<T>::operator() (int r, int c)
{
	int i = getConnection(r, c);
	UG_ASSERT(i!=-1, "connection ("<< r << ", " << c << ") not there");
	return pRowStart[r][i].dValue;
}

////////////////////////////////////////////////////////////////////////////////////



// add Submatrix
//--------------------
//! function to add submatrices ( submatrix )
template<typename T>
template<typename M>
void SparseMatrix<T>::add(const M &mat, size_t *rows, size_t *cols)
{
	connection *c = new connection[mat.num_cols()];
	size_t nc;
	for(size_t i=0; i < mat.num_rows(); i++)
	{
		nc = 0;
		for(size_t j=0; j < mat.num_cols(); j++)
		{
			if(mat(i,j) != 0.0)
			{
				c[nc].iIndex = cols[j];
				c[nc].dValue = mat(i, j);
				nc++;
			}
		}
		UG_ASSERT(nc < mat.num_cols(), "???");

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
	connection *c = new connection[mat.num_cols()];
	size_t nc;
	for(size_t i=0; i < mat.num_rows(); i++)
	{
		nc = 0;
		for(size_t j=0; j < mat.num_cols(); j++)
		{
			if(mat(i,j) != 0.0)
			{
				c[nc].iIndex = cols[j];
				c[nc].dValue = mat(i, j);
				nc++;
			}
		}
		UG_ASSERT(nc < mat.num_cols(), "???");
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
	vector<sortStruct<size_t> > sortedCols(mat.num_cols());

	for(size_t i=0; i<mat.num_cols(); i++)
	{
		sortedCols[i].index = i;
		sortedCols[i].sortValue = cols[i];
	}
	sort(sortedCols.begin(), sortedCols.end());

	for(size_t i=0; i < mat.num_rows(); i++)
	{
		size_t iindex_global = rows[i];
		size_t iindex_local = i;

		cRowIterator conn = beginRow(iindex_global);
		// diagonal

		size_t j = 0;
		while(j < mat.num_cols() && !conn.isEnd())
		{
			size_t jindex_global = sortedCols[j].sortValue;
			size_t jindex_local = sortedCols[j].index;

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
	ASSERT1(mat.num_cols() == cols.size() && mat.num_rows() == rows.size());
	add(mat, &rows[0], &cols[0]);
}

template<typename T>
template<typename M>
void SparseMatrix<T>::set(const M &mat, vector<size_t> &rows, vector<size_t> &cols)
{
	ASSERT1(mat.num_cols() == cols.size() && mat.num_rows() == rows.size());
	set(mat, &rows[0], &cols[0]);
}

template<typename T>
template<typename M>
void SparseMatrix<T>::get(M &mat, vector<size_t> &rows, vector<size_t> &cols) const
{
	ASSERT1(mat.num_cols() == cols.size() && mat.num_rows() == rows.size());
	get(mat, &rows[0], &cols[0]);
}

//////////////////////////////////////////////////////////////////////////

#ifndef FLEXAMG
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
#endif FLEXAMG
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



} // namespace ug
