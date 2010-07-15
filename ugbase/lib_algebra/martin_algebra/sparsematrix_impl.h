/*
 *  sparsematrix_impl.hpp
 *  flexamg
 *
 *  Created by Martin Rupp on 04.11.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 */
#ifndef __H__UG__MARTIN_ALGEBRA__SPARSEMATRIX_IMPL__
#define __H__UG__MARTIN_ALGEBRA__SPARSEMATRIX_IMPL__

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

//======================================================================================================
// construction etc.

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

// resize
//---------------------
//! used to resize the SparseMatrix
//! @param newRows new nr of rows
//! @param newCols new nr of cols
template<typename T>
bool SparseMatrix<T>::resize(size_t newRows, size_t newCols)
{
	if(cols == 0 && rows == 0)
		return create(newRows, newCols);
	if(newCols < cols)
	{
		definalize();
		// remove all connections A(r, c) with c >= newCols
		for(int r=0; r < rows; r++)
		{
			size_t nr;
			if(get_connection_nr(r, newCols, nr, GREATER_EQUAL))
				pRowEnd[r] = pRowStart[r]+nr;
		}
	}

	if(newRows < rows)
	{
		// safe delete rows with r >= newRows.
		for(int r=newRows; r<rows; r++)
			safe_set_connections(r, NULL);
	}

	if(newRows != rows)
	{
		// reallocate arrays
		connection **pNewRowStart = new connection*[newRows+1];
		memcpy(pNewRowStart, pRowStart, sizeof(connection*)*(newRows+1));
		delete[] pRowStart;
		pRowStart = pNewRowStart;

		if(is_finalized())
			pRowEnd = pRowStart+1;
		else
		{
			connection **pNewRowEnd = new connection*[newRows+1];
			memcpy(pNewRowEnd, pRowEnd, sizeof(connection*)*(newRows+1));
			delete[] pRowEnd;
			pRowEnd = pNewRowEnd;
		}

		int *iNewMaxNrOfConnections = new size_t[newRows];
		memcpy(iNewMaxNrOfConnections, iMaxNrOfConnections, sizeof(size_t)*newRows);
		delete[] iMaxNrOfConnections;
		iMaxNrOfConnections = iNewMaxNrOfConnections;
	}

	return true;
}

template<typename T> bool
SparseMatrix<T>::destroy()
{
	if(is_finalized())
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
				safe_set_connections(i, NULL);
			delete [] pRowStart;
			pRowStart = NULL;
		}
		if(pRowEnd) delete [] pRowStart; pRowStart = NULL;
	}

	if(consmem)	delete [] consmem; consmem = NULL;
	rows = cols = 0;
	return true;
}

// createAsTransposeOf
//-----------------------
//!
//! write in a empty SparseMatrix the transpose SparseMatrix of B.
template<typename T>
void SparseMatrix<T>::create_as_transpose_of(const SparseMatrix &B)
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
			
			UG_ASSERT(k>=0 && k<num_connections(ndx), "k = " << k << ". precalculated nr of Connections " << num_connections(ndx) << " wrong?");
			
			pRowStart[ndx][k].dValue = (*conn).dValue;
			pRowStart[ndx][k].iIndex = i;
			if(bandwidth < i-ndx) bandwidth = i-ndx;
			nr[ndx]++;
		}
	
	// note that all connections are already sorted.
	delete[] nr;
}

//======================================================================================================
// finalizing functions


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// defrag
/** defragmentates the matrix by writing all matrix rows consecutively in memory.
 * Sets pRowEnd = pRowStart+1.
 */
template<typename T>
void SparseMatrix<T>::defrag()
{
	iTotalNrOfConnections=0;
	for(size_t i=0; i<rows; i++)
		iTotalNrOfConnections+= num_connections(i);

	connection *consmemNew = new connection[iTotalNrOfConnections];
	connection *p = consmemNew;
	for(size_t i=0; i<rows; i++)
	{
		size_t nr=num_connections(i);
		for(size_t k=0; k < nr; k++)
			swap(p[k], pRowStart[i][k]);
		safe_set_connections(i, p);
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

////////////////////////////////////////////////////////////////////////////////////////////////////////
// definalize
/** creates own pRowEnd array, so that rows dont have to be consecutive anymore
 */
template<typename T>
void SparseMatrix<T>::definalize()
{
	if(!is_finalized())
		return; // already not finalized.
	connection **pNewRowEnd = new connection*[rows];
	memcpy(pNewRowEnd, pRowEnd, sizeof(connection*)*rows);
	pRowEnd = pNewRowEnd;
}

template<typename T>
bool SparseMatrix<T>::is_finalized() const
{
	return pRowEnd == (pRowStart+1);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
// safe_set_connections
/** "safe" way to set a connection, since when cons[row] is
 * in the big consecutive consmem-array, you mustn'd delete it.
 * @param row row to set
 * @param connection pointer to connection to set
 */
template<typename T>
void SparseMatrix<T>::safe_set_connections(size_t row, connection *mem) const
{
	if(pRowStart[row] != NULL && !in_consmem(row))
		delete[] pRowStart[row];
	pRowStart[row] = mem;
}


template<typename T>
bool SparseMatrix<T>::in_consmem(size_t row) const
{
	return pRowStart[row] >= consmem && pRowEnd[row] < consmem + consmemsize;
}





//======================================================================================================
// general functions


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
		if(is_unconnected(i)) continue;
		for(rowIterator conn = beginRow(i); !conn.isEnd(); ++conn)
		{
			size_t conindex = (*conn).iIndex;
			if(is_unconnected(conindex))
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
	if(num_connections(row) == 0 || !in_consmem(row))
	{
		if(in_consmem(row))
			iFragmentedMem += 1;
		else
			iFragmentedMem += 1-num_connections(row);

		connection *c = new connection[1];
		safe_set_connections(row, c);
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

//======================================================================================================
// submatrix set/get



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
			add_matrix_row(rows[i], c, nc);
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
			set_matrix_row(rows[i], c, nc);
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

////////////////////////////////////////////////////////////////////////

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
			add_matrix_row(I[i][0], c, nc);
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
			set_matrix_row(I[i][0], c, nc);
	}
	delete[] c;
	return true;
}
#endif // FLEXAMG

// getDiag
//-------------
//! get Diagonal A_[i,i] of matrix
//! @param i
//! @return A_{i,i}
template<typename T>
inline const T &SparseMatrix<T>::get_diag(size_t i) const
{
	return operator() (i, i);
}

template<typename T>
inline T &SparseMatrix<T>::get_diag(size_t i)
{
	return operator() (i, i);
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
		UG_ASSERT(0, "not tested yet");
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





//======================================================================================================
// row functions


// is_unconnected
//-----------------
/*!
 @remark since first connection is always i,i, every row has at least 1 connection
 @return true if row i has no connection to indices other than i
 */
template<typename T>
inline bool SparseMatrix<T>::is_unconnected(size_t i) const
{
	UG_ASSERT(i < rows && i >= 0, *this << ": " << i << " out of bounds.");
	if(pRowStart[i] == NULL)
		return true;
	int nr=num_connections(i);
	if(nr == 0 || (pRowStart[i][0].iIndex == i && nr == 1))
		return true;
	else
		return false;
}

template<typename T>
inline size_t SparseMatrix<T>::num_connections(size_t row) const
{
	return pRowEnd[row] - pRowStart[row];
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


////////////////////////////////////////////////////////////////////////////////////////////////////
// set_matrix_row
/** set a row of the matrix. all previous content in this row is destroyed (@sa add_matrix_row).
 * @param row index of the row to set
 * @param c pointer to a array of sorted connections of size nr
 * @param nr number of connections in c
 * @remark connections have to be sorted
 */
template<typename T>
void SparseMatrix<T>::set_matrix_row(size_t row, connection *c, size_t nr)
{
	definalize();

	if(in_consmem(row))
		iFragmentedMem += nr;
	else
		iFragmentedMem += nr - num_connections(row);

	connection *n = new connection[nr];
	for(size_t i=0; i<nr; i++)
		n[i] = c[i];

	sort(n, n+nr);

	for(size_t i=0; i<nr; i++)
		UG_ASSERT(n[i].iIndex >= 0 && n[i].iIndex < num_cols(), "Matrix is " << num_rows() << "x" << num_cols() << ", row " << row << " cannot have connection " << n[i] << ".");


	bandwidth = max(bandwidth, abs(n[0].iIndex, row));
	bandwidth = max(bandwidth, abs(n[nr-1].iIndex, row));
	

	definalize();
	safe_set_connections(row, n);
	pRowEnd[row] = pRowStart[row]+nr;

	iTotalNrOfConnections += nr - num_connections(row);


	iMaxNrOfConnections[row] = nr;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// add_matrix_row
/** add a row to a matrix row.
 * @param row index of the row to set
 * @param c pointer to a array of sorted connections of size nr
 * @param nr number of connections in c
 * @remark if we get new connections, matrix is definalized.
 * @remark connections have to be sorted
 */
template<typename T>
void SparseMatrix<T>::add_matrix_row(size_t row, connection *c, size_t nr)
{	
	if(nr <= 0) return;

	//cout << "row: " << row << " nr: " << nr << endl;
	connection *old = pRowStart[row];

	if(old == NULL)
	{		
		set_matrix_row(row, c, nr);
		return;
	}
	
	UG_ASSERT(num_connections(row) != 0, "cons[row] != NULL but get_nr_of_connection(row) == 0 ???");

	// the matrix row is not empty and we are adding more than one connection

	size_t oldNrOfConnections = num_connections(row);

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

	if(in_consmem(row))
		iFragmentedMem += iNewSize;
	else
		iFragmentedMem += iNewSize - num_connections(row);

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
	safe_set_connections(row, n);
	pRowEnd[row] = pRowStart[row]+iNewSize;

	iTotalNrOfConnections += iNewSize - oldNrOfConnections;
	iMaxNrOfConnections[row] = iNewSize;
}


///////////////////////////////////////////////////////////////////////////////


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
// <=, <, =, >=, =



//======================================================================================================
// connectivity functions.

////////////////////////////////////////////////////////////////////////////////////////////////////
// get_connection
/** returns a rowIterator to the connection A(r,c), creates connection if necessary.
 * @param r index of the row
 * @param c index of the column
 * @remark creates connection if necessary.
 */
template<typename T>
typename SparseMatrix<T>::rowIterator SparseMatrix<T>::get_connection(size_t r, size_t c)
{
	size_t nr;
	bool bFound = get_connection_nr(r, c, nr, GREATER_EQUAL);

	if(!bFound || pRowStart[r][nr].iIndex != c)
	{
		int numConnections = num_connections(r);
		if(bFound == false)
			nr = numConnections;
		connection *con = new connection[numConnections+1];
		if(pRowStart[r]) memcpy(con, pRowStart[r], nr*sizeof(connection));
		con[nr].iIndex = c;
		con[nr].dValue = 0.0;
		if(pRowStart[r]) memcpy(con+nr+1, pRowStart[r], (numConnections-nr)*sizeof(connection));
		definalize();
		safe_set_connections(r, con);
		pRowEnd[r] = pRowStart[r]+nr+1;
	}

	rowIterator it=beginRow(r);
	it += nr;
	return it;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// get_connection
/** returns a cRowIterator to the connection A(r,c) if connection already there
 * @param r index of the row
 * @param c index of the column
 * @param bFound returns true if connection found, otherwise false.
 * @remark does not create connections
 */
template<typename T>
typename SparseMatrix<T>::cRowIterator SparseMatrix<T>::get_connection(size_t r, size_t c, bool &bFound) const
{
	cRowIterator it=beginRow(r);

	size_t nr;
	bFound = get_connection_nr(r, c, nr);
	if(bFound == false)
		it.p = it.pEnd;
	else
		it.p += nr;

	return it;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// get_connection
/** returns a rowIterator to the connection A(r,c) if connection already there
 * @param r index of the row
 * @param c index of the column
 * @param bFound returns true if connection found, otherwise false.
 * @remark does not create connections
 */
template<typename T>
typename SparseMatrix<T>::rowIterator SparseMatrix<T>::get_connection(size_t r, size_t c, bool &bFound)
{
	rowIterator it=beginRow(r);

	size_t nr;
	bFound = get_connection_nr(r, c, nr);
	if(bFound == false)
		it.p = it.pEnd;
	else
		it.p += nr;

	return it;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// get_connection_nr_templ
/** binary searches connection A(r,c).
 * @param r index of the row
 * @param c index of the column
 * @param nr returns nr in pRowStart[r], so that pRowStart[r][nr].iIndex = c or >=/>/</<= @sa get_connection_nr
 * @return true if connection was exactly found in = mode, otherwise always true.
 * @remarks >=/>/=/<=/< depends on template parameters
 */
template<typename T>
template<typename SparseMatrix<T>::get_connection_nr_flag flag>
bool SparseMatrix<T>::get_connection_nr_templ(size_t r, size_t c, size_t &nr) const
{
	nr = 0;
	connection *con = pRowStart[r];
	UG_ASSERT(con != NULL, "");

	if(num_connections(r) == 0) return false;
	if(num_connections(r) == 1 && flag == EQUAL)
		return con[0].iIndex == c;

	size_t left=0, right = num_connections(r)-1;

	while(left < right)
	{
		size_t mid = (left+right)/2;
		if (con[mid].iIndex < c)
			left = mid+1;
		else if (con[mid].iIndex > c)
			right = mid-1;
		else
		{
			if(flag == LESS_EQUAL || flag == EQUAL || flag == GREATER_EQUAL)
			{
				nr = mid;
				return true;
			}
			left = right = mid;
			break;
		}
	}

	if((flag == LESS_EQUAL || flag == EQUAL || flag == GREATER_EQUAL) && left == right && con[left].iIndex == c)
	{
		nr = left;
		return true;
	}
	else if(flag == LESS || flag == LESS_EQUAL)
	{
		if(left == 0)
		{
			if(con[0].iIndex < c)
			{
				nr = 0;
				return true;
			}
			return false;
		}

		nr = left-1;
		return true;
	}
	else if(flag == GREATER || flag == GREATER_EQUAL)
	{
		if(right+1 >= num_connections(r))
		{
			if(con[num_connections(r)-1].iIndex > c)
			{
				nr = num_connections(r)-1;
				return true;
			}
			return false;
		}

		nr = right+1;
		return true;
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// get_connection_nr
/** searches connections A(r,c)
 * @param r index of the row
 * @param c index of the column
 * @param nr returns nr in pRowStart[r], so that pRowStart[r][nr].iIndex = c or >=/>/</<= depending on flag
 * @param flag EQUAL, LESS_EQUAL, LESS, GREATER, or GREATER_EQUAL
 * @return true if connection was exactly found in = mode, otherwise always false if no connection at all.
 * @remarks flag=EQUAL is standard parameter
 */
template<typename T>
bool SparseMatrix<T>::get_connection_nr(size_t r, size_t c, size_t &nr, get_connection_nr_flag flag) const
{
	if(pRowStart[r] == NULL) return false;
	switch(flag)
	{
	case LESS: return get_connection_nr_templ<LESS>(r, c, nr);
	case LESS_EQUAL: return get_connection_nr_templ<LESS_EQUAL>(r, c, nr);
	case EQUAL: return get_connection_nr_templ<EQUAL>(r, c, nr);
	case GREATER_EQUAL: return get_connection_nr_templ<GREATER_EQUAL>(r, c, nr);
	case GREATER: return get_connection_nr_templ<GREATER>(r, c, nr);
	default: UG_ASSERT(0, "non-existing flag option"); return false;
	}
}


template<typename T>
const typename SparseMatrix<T>::entry_type &SparseMatrix<T>::operator() (size_t r, size_t c, bool &bFound) const
{
	size_t nr=0;
	bFound = get_connection_nr(r, c, nr);
	if(!bFound)
	{
		bFound = false;
		static T t; t=0.0;
		return t;
	}
	else
		return pRowStart[r][nr].dValue;
}

template<typename T>
typename SparseMatrix<T>::entry_type &SparseMatrix<T>::operator() (size_t r, size_t c, bool &bFound)
{
	size_t nr=0;
	bFound = get_connection_nr(r, c, nr);
	if(!bFound)
	{
		bFound = false;
		static T t;	t = 0.0;
		return t;
	}
	else
		return pRowStart[r][nr].dValue;
}

template<typename T>
const typename SparseMatrix<T>::entry_type &SparseMatrix<T>::operator() (size_t r, size_t c) const
{
	size_t nr=0;
	bool bConnectionFound = get_connection_nr(r, c, nr);
	if(bConnectionFound != true)
		UG_ASSERT(bConnectionFound == true, "connection A(" << r << "," << c << ") not there and A const. Use operator(r, c, bFound) to catch.");
	return pRowStart[r][nr].dValue;
}

template<typename T>
typename SparseMatrix<T>::entry_type &SparseMatrix<T>::operator() (size_t r, size_t c)
{
	return (*get_connection(r, c)).dValue;
}



} // namespace ug

#endif

