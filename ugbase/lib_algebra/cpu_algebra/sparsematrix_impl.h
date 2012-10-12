/*
 *  sparsematrix_impl.hpp
 *
 *  Created by Martin Rupp on 04.11.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 */
#ifndef __H__UG__CPU_ALGEBRA__SPARSEMATRIX_IMPL__
#define __H__UG__CPU_ALGEBRA__SPARSEMATRIX_IMPL__

// creation etc

#include <fstream>
#include <cstring>

#include "algebra_misc.h"
#include "local_helper.h"
#include "lib_algebra/common/operations_vec.h"
#include "common/profiler/profiler.h"
#include <vector>
#include <algorithm>

template<typename T> T abs(const T &a, const T &b)
{
	if(a > b) return a-b; else return b-a;
}

namespace ug{

//======================================================================================================
// construction etc.

template<typename T>
SparseMatrix<T>::SparseMatrix() : rows(0), cols(0), pRowStart(NULL), pRowEnd(NULL),
	iTotalNrOfConnections(0), bandwidth(0), estimatedRowSize(0), iMaxNrOfConnections(0), consmem(NULL),
	 consmemsize(0), iFragmentedMem(0), m_bIgnoreZeroes(false)
{
	FORCE_CREATION { print(); p(); pr(0); }
}


template<typename T>
SparseMatrix<T>::~SparseMatrix()
{
	destroy();
}



template<typename T>
bool SparseMatrix<T>::create(size_t _rows, size_t _cols)
{
	UG_ASSERT(rows == 0 && cols == 0, *this << " not empty.");

	rows = _rows;
	cols = _cols;

	pRowStart = new connection*[rows+1];
	UG_ASSERT(pRowStart != NULL, "out of memory, no more space for " << sizeof(connection*)*(rows+1));
	memset(pRowStart, 0, sizeof(connection*)*(rows+1));

	pRowEnd = new connection*[rows+1];
	UG_ASSERT(pRowEnd != NULL, "out of memory, no more space for " << sizeof(connection*)*(rows+1));
	memset(pRowEnd, 0, sizeof(connection*)*(rows+1));

	if(rows == 0)
		iMaxNrOfConnections = NULL;
	else {
		iMaxNrOfConnections = new size_t[rows];
		UG_ASSERT(iMaxNrOfConnections != NULL, "out of memory, no more space for " << sizeof(size_t)*rows);
		memset(iMaxNrOfConnections, 0, sizeof(size_t)*rows);
	}

	iTotalNrOfConnections = 0;
	bandwidth = 0;
	return true;
}


template<typename T>
bool SparseMatrix<T>::resize(size_t newRows, size_t newCols)
{
	PROFILE_BEGIN_GROUP(SparseMatrix_resize, "algebra SparseMatrix");

	if(cols == 0 && rows == 0)
	{
		// 	\todo: verify: A matrix could not be created, but is resized
		if(is_finalized())
		{
			if(pRowStart) delete [] pRowStart;
			pRowStart = NULL;
			pRowEnd = NULL;
		}
		else
		{
			if(pRowStart) delete[] pRowStart;
			if(pRowEnd) delete[] pRowEnd;
		}
		return create(newRows, newCols);
	}
	if(newCols < cols)
	{
		definalize();
		// remove all connections A(r, c) with c >= newCols
		for(size_t r=0; r < rows; r++)
		{
			size_t nr;
			if(get_connection_nr(r, newCols, nr, GREATER_EQUAL))
			{
				iTotalNrOfConnections -= num_connections(r)-nr;
				pRowEnd[r] = pRowStart[r]+nr;
			}
		}
		cols = newCols;
	}
	else
		cols = newCols;

	if(newRows < rows)
	{
		// safe delete rows with r >= newRows.
		for(size_t r=newRows; r<rows; r++)
		{
			iTotalNrOfConnections -= num_connections(r);
			safe_set_connections(r, NULL);
		}
	}

	if(newRows != rows)
	{
		bool bWasFinalized = is_finalized();
		// reallocate arrays
		connection **pNewRowStart = new connection*[newRows+1];
		UG_ASSERT(pNewRowStart != NULL, "out of memory, no more space for " << sizeof(connection*)*(newRows+1));
		memcpy(pNewRowStart, pRowStart, sizeof(connection*)*std::min(rows+1, newRows+1));
		delete[] pRowStart;
		pRowStart = pNewRowStart;

		if(bWasFinalized)
		{
			pRowEnd = pRowStart+1;
			for(size_t i=rows+1; i<newRows+1; i++)
				pRowStart[i] = pRowStart[rows];
		}
		else
		{
			connection **pNewRowEnd = new connection*[newRows+1];
			UG_ASSERT(pNewRowEnd != NULL, "out of memory, no more space for " << sizeof(connection*)*(newRows+1));
			memcpy(pNewRowEnd, pRowEnd, sizeof(connection*) * std::min(rows, newRows));
			delete[] pRowEnd;
			pRowEnd = pNewRowEnd;

			for(size_t i=rows; i<newRows; i++)
			{
				pRowStart[i] = NULL;
				pRowEnd[i] = NULL;
			}
			pRowStart[newRows] = NULL;
		}

		if(newRows > 0) {
			size_t *iNewMaxNrOfConnections = new size_t[newRows];

			UG_ASSERT(iNewMaxNrOfConnections != NULL, "out of memory, no more space for " << sizeof(size_t)*newRows);
			memcpy(iNewMaxNrOfConnections, iMaxNrOfConnections, sizeof(size_t) * std::min(rows, newRows));
			delete[] iMaxNrOfConnections;
			iMaxNrOfConnections = iNewMaxNrOfConnections;

			for(size_t i=rows; i<newRows; i++)
				iMaxNrOfConnections[i] = 0;
		}

		rows = newRows;
	}

	return true;
}

template<typename T> bool
SparseMatrix<T>::destroy()
{
	PROFILE_BEGIN_GROUP(SparseMatrix_destroy, "algebra SparseMatrix");
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
		if(pRowEnd) delete [] pRowEnd; pRowEnd = NULL;
	}

	if(consmem)	delete [] consmem; consmem = NULL;
	if(iMaxNrOfConnections) delete[] iMaxNrOfConnections; iMaxNrOfConnections = NULL;
	rows = cols = 0;
	return true;
}

template<typename T>
bool SparseMatrix<T>::set_as_transpose_of(const SparseMatrix<T> &B, double scale)
{
	PROFILE_BEGIN_GROUP(SparseMatrix_set_as_transpose_of, "algebra SparseMatrix");
	//destroy();
	UG_ASSERT(rows == 0 && cols == 0, *this << " not empty.");

	rows = B.cols;
	cols = B.rows;

	pRowStart = new connection*[rows+1];
	UG_ASSERT(pRowStart != NULL, "out of memory, no more space for " << sizeof(connection*)*(rows+1));

	pRowEnd= pRowStart+1;
	iMaxNrOfConnections = new size_t[rows];
	UG_ASSERT(iMaxNrOfConnections != NULL, "out of memory, no more space for " << sizeof(size_t)*rows);

	iTotalNrOfConnections = 0;
	bandwidth = 0;

	//-----------------------------------

	size_t *nr = new size_t[rows];
	UG_ASSERT(nr != NULL, "out of memory, no more space for " << sizeof(size_t)*rows);
	memset(nr, 0, sizeof(size_t)*rows);

	// get length of each row
	for(size_t j=0; j < B.rows; j++)
		for(const_row_iterator conn = B.begin_row(j); conn != B.end_row(j); ++conn)
			if(conn.value() == 0) continue;
			else nr[conn.index()]++;


	size_t newTotal = 0;
	for(size_t i=0; i < rows; i++)
		newTotal += nr[i];


	// allocate one big connection memory
	consmem = new connection[newTotal];
	UG_ASSERT(consmem != NULL, "out of memory, no more space for " << sizeof(connection)*newTotal);
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

	if(scale == 1.0)
	{
		for(size_t i=0; i < B.rows; i++)
			for(const_row_iterator conn = B.begin_row(i); conn != B.end_row(i); ++conn)
			{
				if(conn.value() == 0) continue;
				size_t ndx = conn.index();
				UG_ASSERT(ndx >= 0 && ndx < rows, "connection " << (*conn) << " of " << B << ", row " << i << " out of range 0.." << rows);

				size_t k = nr[ndx];

				UG_ASSERT(k>=0 && k<num_connections(ndx), "k = " << k << ". precalculated nr of Connections " << num_connections(ndx) << " wrong?");

				pRowStart[ndx][k].dValue = conn.value();
				pRowStart[ndx][k].iIndex = i;
				if(bandwidth < i-ndx) bandwidth = i-ndx;
				nr[ndx]++;
			}
	}
	else
	{
		for(size_t i=0; i < B.rows; i++)
			for(const_row_iterator conn = B.begin_row(i); conn != B.end_row(i); ++conn)
			{
				if(conn.value() == 0) continue;
				size_t ndx = conn.index();
				UG_ASSERT(ndx >= 0 && ndx < rows, "connection " << (*conn) << " of " << B << ", row " << i << " out of range 0.." << rows);

				size_t k = nr[ndx];

				UG_ASSERT(k>=0 && k<num_connections(ndx), "k = " << k << ". precalculated nr of Connections " << num_connections(ndx) << " wrong?");

				pRowStart[ndx][k].dValue = conn.value();
				pRowStart[ndx][k].dValue *= scale;
				pRowStart[ndx][k].iIndex = i;
				if(bandwidth < i-ndx) bandwidth = i-ndx;
				nr[ndx]++;
			}
	}


	// note that all connections are already sorted.
	delete[] nr;
	return true;
}

//======================================================================================================
// finalizing functions


/////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
void SparseMatrix<T>::defragment()
{
	if(is_finalized()) return;
	PROFILE_BEGIN_GROUP(SparseMatrix_defragment, "algebra SparseMatrix");

	iTotalNrOfConnections=0;
	for(size_t i=0; i<rows; i++)
		iTotalNrOfConnections += num_connections(i);

	connection *consmemNew = new connection[iTotalNrOfConnections+1];
	UG_ASSERT(consmemNew != NULL, "out of memory, no more space for " << sizeof(connection)*iTotalNrOfConnections);
	connection *p = consmemNew;
	for(size_t i=0; i<rows; i++)
	{
		size_t nr=num_connections(i);
		for(size_t k=0; k < nr; k++)
			std::swap(p[k], pRowStart[i][k]);
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

template<typename T>
void SparseMatrix<T>::definalize()
{
	if(!is_finalized())
		return; // already not finalized.
	connection **pNewRowEnd = new connection*[rows];
	UG_ASSERT(pNewRowEnd != NULL, "out of memory, no more space for " << sizeof(connection*)*rows);
	memcpy(pNewRowEnd, pRowEnd, sizeof(connection*)*rows);
	pRowEnd = pNewRowEnd;
	// no need to delete old pRowEnd, since it was pRowEnd = pRowStart+1.
}

template<typename T>
bool SparseMatrix<T>::is_finalized() const
{
	return pRowEnd == (pRowStart+1);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
void SparseMatrix<T>::safe_set_connections(size_t row, connection *mem) const
{
	UG_ASSERT(row < num_rows(), "cannot access row " << row << " of " << *this);

	if(pRowStart[row] != NULL && !in_consmem(row))
		delete[] pRowStart[row];
	pRowStart[row] = mem;
}


template<typename T>
bool SparseMatrix<T>::in_consmem(size_t row) const
{
	UG_ASSERT(row < num_rows(), "cannot access row " << row << " of " << *this);
	return pRowStart[row] >= consmem && pRowEnd[row] <= consmem + consmemsize;
}


//======================================================================================================
// general functions


// calculate dest = alpha1*v1 + beta1*A*w1 (A = this matrix)
template<typename T>
template<typename vector_t>
bool SparseMatrix<T>::axpy(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const vector_t &w1) const
{
	PROFILE_BEGIN_GROUP(SparseMatrix_axpy, "algebra SparseMatrix");
//	UG_ASSERT(cols == x.size(), "x: " << x << " has wrong length (should be " << cols << "). A: " << *this);
//	UG_ASSERT(rows == res.size(), "res: " << x << " has wrong length (should be " << rows << "). A: " << *this);


	if(alpha1 == 0.0)
	{
		for(size_t i=0; i < rows; i++)
		{
			const_row_iterator conn = begin_row(i);
			if(conn  == end_row(i)) continue;
			MatMult(dest[i], beta1, conn.value(), w1[conn.index()]);
			for(++conn; conn != end_row(i); ++conn)
				// res[i] += conn.value() * x[conn.index()];
				MatMultAdd(dest[i], 1.0, dest[i], beta1, conn.value(), w1[conn.index()]);
		}
	}
	else if(&dest == &v1)
	{
		if(alpha1 != 1.0)
			for(size_t i=0; i < rows; i++)
			{
				dest[i] *= alpha1;
				mat_mult_add_row(i, dest[i], beta1, w1);
			}
		else
			for(size_t i=0; i < rows; i++)
				mat_mult_add_row(i, dest[i], beta1, w1);

	}
	else
	{
		for(size_t i=0; i < rows; i++)
		{
			VecScaleAssign(dest[i], alpha1, v1[i]);
			mat_mult_add_row(i, dest[i], beta1, w1);
		}
	}

	return true;
}

// calculate dest = alpha1*v1 + beta1*A^T*w1 (A = this matrix)
template<typename T>
template<typename vector_t>
bool SparseMatrix<T>::axpy_transposed(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const vector_t &w1) const
{
	PROFILE_BEGIN_GROUP(SparseMatrix_axpy_transposed, "algebra SparseMatrix");
//	UG_ASSERT(rows == x.size(), "x: " << x << " has wrong length (should be " << rows << "). A: " << *this);
//	UG_ASSERT(cols == res.size(), "res: " << x << " has wrong length (should be " << cols << "). A: " << *this);

	if(&dest == &v1)
	{
		if(alpha1 != 1.0)
		{
			if(alpha1 == 0.0)
				dest.set(0.0);
			else if(alpha1 != 1.0)
				dest *= alpha1;
		}
	}
	else if(alpha1 != 0.0)
		VecScaleAssign(dest, alpha1, v1);
	else
		dest.set(0.0);

	for(size_t i=0; i<rows; i++)
	{
		for(const_row_iterator conn = begin_row(i); conn != end_row(i); ++conn)
		{
			if(conn.value() != 0.0)
				// dest[conn.index()] += beta1 * conn.value() * w1[i];
				MatMultTransposedAdd(dest[conn.index()], 1.0, dest[conn.index()], beta1, conn.value(), w1[i]);
		}
	}
	return true;
}

template<typename T>
template<typename vector_t>
inline void SparseMatrix<T>::mat_mult_add_row(size_t row, typename vector_t::value_type &dest, double alpha, const vector_t &v) const
{
	for(const_row_iterator conn = begin_row(row); conn != end_row(row); ++conn)
		MatMultAdd(dest, 1.0, dest, alpha, conn.value(), v[conn.index()]);
}


//======================================================================================================
// submatrix set/get

template<typename T>
template<typename M>
void SparseMatrix<T>::add(const M &mat)
{
	connection *c = new connection[mat.num_cols()];
	UG_ASSERT(c != NULL, "out of memory, no more space for " << sizeof(connection)*mat.num_cols());

	size_t nc;
	for(size_t i=0; i < mat.num_rows(); i++)
	{
		nc = 0;
		for(size_t j=0; j < mat.num_cols(); j++)
		{
			if(!m_bIgnoreZeroes || mat(i,j) != 0.0)
			{
				c[nc].iIndex = mat.col_index(j);
				c[nc].dValue = mat(i, j);
				nc++;
			}
		}
		UG_ASSERT(nc <= mat.num_cols(), "More column entries read, than local matrix contains.");

		if(nc > 0)
			add_matrix_row(mat.row_index(i), c, nc);
	}
	delete[] c;
}


template<typename T>
template<typename M>
void SparseMatrix<T>::set(const M &mat)
{
	connection *c = new connection[mat.num_cols()];
	UG_ASSERT(c != NULL, "out of memory, no more space for " << sizeof(connection)*mat.num_cols());

	size_t nc;
	for(size_t i=0; i < mat.num_rows(); i++)
	{
		nc = 0;
		for(size_t j=0; j < mat.num_cols(); j++)
		{
			if(!m_bIgnoreZeroes || mat(i,j) != 0.0)
			{
				c[nc].iIndex = mat.col_index(j);
				c[nc].dValue = mat(i, j);
				nc++;
			}
		}
		UG_ASSERT(nc <= mat.num_cols(), "More column entries read, than local matrix contains.");

		if(nc > 0)
			set_matrix_row(mat.row_index(i), c, nc);
	}
	delete[] c;
}

template<typename T>
template<typename M>
void SparseMatrix<T>::get(M &mat) const
{
	std::vector<sortStruct<size_t> > sortedCols(mat.num_cols());

	for(size_t i=0; i<mat.num_cols(); i++)
	{
		sortedCols[i].index = i;
		sortedCols[i].sortValue = mat.col_index(i);
	}
	std::sort(sortedCols.begin(), sortedCols.end());

	for(size_t i=0; i < mat.num_rows(); i++)
	{
		size_t iindex_global = mat.row_index(i);
		size_t iindex_local = i;

		const_row_iterator conn = begin_row(iindex_global);
		// diagonal

		size_t j = 0;
		while(j < mat.num_cols() && conn != end_row(iindex_global))
		{
			size_t jindex_global = sortedCols[j].sortValue;
			size_t jindex_local = sortedCols[j].index;

			size_t cindex_global = conn.index();

			if(cindex_global < jindex_global)
				++conn;
			else if(cindex_global > jindex_global)
			{
				mat(iindex_local, jindex_local) = 0.0;
				++j;
			}
			else
			{
				mat(iindex_local, jindex_local) = conn.value();
				++conn; ++j;
			}
		}
	}
}
/*
template<typename T>
template<typename M>
void SparseMatrix<T>::add(const M &mat, size_t *rows, size_t *cols)
{
	add(c_localMatrix_from_mat_and_array<M> (mat, rows, cols));
}

template<typename T>
template<typename M>
void SparseMatrix<T>::set(const M &mat, size_t *rows, size_t *cols)
{
	set(c_localMatrix_from_mat_and_array<M> (mat, rows, cols));
}

template<typename T>
template<typename M>
void SparseMatrix<T>::get(M &mat, size_t *rows, size_t *cols) const
{
	localMatrix_from_mat_and_array<M> g(mat, rows, cols);
	get(g);
}*/



template<typename T>
bool SparseMatrix<T>::set(double a)
{
	if(a == 0.0)
	{
		UG_ASSERT(pRowStart != NULL, "");
		for(size_t row=0; row<rows; row++)
			for(row_iterator it = begin_row(row); it != end_row(row); ++it)
				it.value() = 0.0;
	}
	else
	{
		UG_ASSERT(0, "not tested yet");
		for(size_t row=0; row<rows; row++)
			for(row_iterator it = begin_row(row); it != end_row(row); ++it)
			{
				if(it.index() == row)
					it.value() = a;
				else
					it.value() = 0.0;
			}
	}
	return true;
}





//======================================================================================================
// row functions


template<typename T>
inline bool SparseMatrix<T>::is_isolated(size_t i) const
{
	UG_ASSERT(i < rows && i >= 0, *this << ": " << i << " out of bounds.");
	if(pRowStart[i] == NULL)
		return true;

	for(const_row_iterator it = begin_row(i); it != end_row(i); ++it)
		if(it.index() != i && it.value() != 0.0)
			return false;
	return true;
}

template<typename T>
inline size_t SparseMatrix<T>::num_connections(size_t row) const
{
	return pRowEnd[row] - pRowStart[row];
}


////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
void SparseMatrix<T>::set_matrix_row(size_t row, connection *c, size_t nr)
{
	if(nr==0)
	{
		safe_set_connections(row, NULL);
		pRowEnd[row] = pRowStart[row];
		return;
	}

	definalize();
	size_t oldNumConnections = num_connections(row);
	if(in_consmem(row))
		iFragmentedMem += nr;
	else
		iFragmentedMem += nr - oldNumConnections;

	std::sort(c, c+nr);
	size_t nrUnique=1;
	for(size_t i=1; i<nr; i++)
	{
		if(c[i].iIndex != c[i-1].iIndex)
			nrUnique++;
	}

	connection *n = new connection[nrUnique];
	UG_ASSERT(n != NULL, "out of memory, no more space for " << sizeof(connection)*nrUnique);

	size_t i=0;
	for(size_t ic=0; ic<nr;)
	{
		n[i] = c[ic++];
		// check for double connections
		while(c[ic-1].iIndex == c[ic].iIndex && ic < nr)
		{
			n[i].dValue += c[ic].dValue;
			ic++;
		}
		i++;
	}

#ifdef DEBUG
	for(size_t i=0; i<nrUnique; i++)
		UG_ASSERT(n[i].iIndex >= 0 && n[i].iIndex < num_cols(), "Matrix is " << num_rows() << "x" << num_cols() << ", row " << row << " cannot have connection " << n[i] << ".");
#endif

	// calc bandwidth
	bandwidth = std::max(bandwidth, abs(n[0].iIndex, row));
	bandwidth = std::max(bandwidth, abs(n[nrUnique-1].iIndex, row));

	safe_set_connections(row, n);
	pRowEnd[row] = pRowStart[row]+nrUnique;

	iTotalNrOfConnections += nrUnique - oldNumConnections;
	iMaxNrOfConnections[row] = nrUnique;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// add_matrix_row
template<typename T>
void SparseMatrix<T>::add_matrix_row(size_t row, connection *c, size_t nr)
{
	if(nr <= 0) return;

	//cout << "row: " << row << " nr: " << nr << endl;
	connection *old = pRowStart[row];

	if(old == NULL || num_connections(row) == 0)
	{
		set_matrix_row(row, c, nr);
		return;
	}

	// the matrix row is not empty and we are adding at least one connection
	size_t oldNrOfConnections = num_connections(row);

	// sort the connections: index1 < index2 < ...
	std::sort(c, c+nr);

	size_t ic=0, iold=0, skipped=0;
	// - add connections which are there
	// - count not existing connections (=skipped)
	while(ic < nr && iold < oldNrOfConnections)
	{
		if(c[ic].iIndex < old[iold].iIndex)
		{
			skipped++;
			ic++;
			// check for double connections
			while(c[ic-1].iIndex == c[ic].iIndex && ic < nr) ic++;
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
	UG_ASSERT(n != NULL, "out of memory, no more space for " << sizeof(connection)*iNewSize);

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

			n[i] = c[ic++];
			// check for double connections
			while(c[ic-1].iIndex == c[ic].iIndex && ic < nr)
			{
				n[i].dValue += c[ic].dValue;
				ic++;
			}
			i++;
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
	while(ic<nr)
	{
		n[i] = c[ic++];
		// check for double connections
		while(c[ic-1].iIndex == c[ic].iIndex && ic < nr)
		{
			n[i].dValue += c[ic].dValue;
			ic++;
		}
		i++;
	}
	while(iold<oldNrOfConnections) n[i++] = old[iold++];

	UG_ASSERT(i == iNewSize, "row: " << row << " i: " << i << " iNewSize: " << iNewSize);

	// set connections
	safe_set_connections(row, n);
	pRowEnd[row] = pRowStart[row]+iNewSize;

	iTotalNrOfConnections += iNewSize - oldNrOfConnections;
	iMaxNrOfConnections[row] = iNewSize;
}


///////////////////////////////////////////////////////////////////////////////


//======================================================================================================
// connectivity functions.

////////////////////////////////////////////////////////////////////////////////////////////////////
// get_connection
template<typename T>
typename SparseMatrix<T>::row_iterator SparseMatrix<T>::get_connection(size_t r, size_t c)
{
	UG_ASSERT(r < num_rows() && c < num_cols(),	"cannot access (" << r << ", " << c << ") of " << *this);

	size_t nr;
	bool bFound = get_connection_nr(r, c, nr, GREATER_EQUAL);

	if(!bFound || pRowStart[r][nr].iIndex != c)
	{
		int numConnections = num_connections(r);
		if(bFound == false)
			nr = numConnections;
		connection *con = new connection[numConnections+1];
		UG_ASSERT(con != NULL, "out of memory, no more space for " << sizeof(connection)*(numConnections+1));

		if(pRowStart[r]) memcpy(con, pRowStart[r], nr*sizeof(connection));
		con[nr].iIndex = c;
		con[nr].dValue = 0.0;
		if(pRowStart[r]) memcpy(con+nr+1, pRowStart[r]+nr, (numConnections-nr)*sizeof(connection));
		definalize();
		safe_set_connections(r, con);
		pRowEnd[r] = pRowStart[r]+numConnections+1;
		iTotalNrOfConnections++;
	}

	row_iterator it = begin_row(r);
	it += nr;
	return it;
}

template<typename T>
typename SparseMatrix<T>::const_row_iterator SparseMatrix<T>::get_connection(size_t r, size_t c) const
{

	UG_ASSERT(r < num_rows() && c < num_cols(),	"cannot access (" << r << ", " << c << ") of " << *this);

	size_t nr=0;
	bool bConnectionFound = get_connection_nr(r, c, nr);
	if(bConnectionFound != true)
		return end_row(r);

	const_row_iterator it = begin_row(r);
	it += nr;
	return it;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
// get_connection
/**
 * returns a const_row_iterator to the connection A(r,c) if connection already there
 * \param r index of the row
 * \param c index of the column
 * \param bFound returns true if connection found, otherwise false.
 * \remark does not create connections
 */
template<typename T>
typename SparseMatrix<T>::const_row_iterator SparseMatrix<T>::get_connection(size_t r, size_t c, bool &bFound) const
{
	UG_ASSERT(r < num_rows() && c < num_cols(),	"cannot access (" << r << ", " << c << ") of " << *this);
	size_t nr;
	bFound = get_connection_nr(r, c, nr);
	if(bFound == false)
		return end_row(r);
	else
	{
		const_row_iterator it=begin_row(r);
		it += nr;
		return it;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// get_connection
/**
 * returns a row_iterator to the connection A(r,c) if connection already there
 * \param r index of the row
 * \param c index of the column
 * \param bFound returns true if connection found, otherwise false.
 * \remark does not create connections
 */
template<typename T>
typename SparseMatrix<T>::row_iterator SparseMatrix<T>::get_connection(size_t r, size_t c, bool &bFound)
{
	UG_ASSERT(r < num_rows() && c < num_cols(),	"cannot access (" << r << ", " << c << ") of " << *this);
	size_t nr;
	bFound = get_connection_nr(r, c, nr);
	if(bFound == false)
		return end_row(r);
	else
	{
		row_iterator it=begin_row(r);
		it += nr;
		return it;
	}
};


////////////////////////////////////////////////////////////////////////////////////////////////////
// get_connection_nr_templ
template<typename T>
template<size_t flag >
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
			{right = mid>0 ? mid-1 : 0;}	// TODO: ATTENTION: I changed this because in case (left-right) == 1 (i.e. mid = left) and left==0, this will lead to huge right: right = (size_t)(-1)
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


		//TODO: PLEASE VERIFY !!! . I change this, because it gave insufficient results.
		//      But i did not rethink all cases with the care I should have taken. Andreas Vogel
		// formerly: nr = right+1;
		if(con[left].iIndex > c)
			nr = left;
		else
			nr = left + 1;
		return true;
	}

	return false;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// get_connection_nr
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
const typename SparseMatrix<T>::value_type &SparseMatrix<T>::operator() (size_t r, size_t c) const
{
	UG_ASSERT(r < num_rows() && c < num_cols(),
				"cannot access (" << r << ", " << c << ") of " << *this);
	size_t nr=0;
	bool bConnectionFound = get_connection_nr(r, c, nr);
	if(bConnectionFound != true)
		UG_ASSERT(bConnectionFound == true, "connection A(" << r << "," << c << ") not there and A const. Use operator(r, c, bFound) to catch.");
	return pRowStart[r][nr].dValue;
}

template<typename T>
typename SparseMatrix<T>::value_type &SparseMatrix<T>::operator() (size_t r, size_t c)
{
	UG_ASSERT(r < num_rows() && c < num_cols(),
				"cannot access (" << r << ", " << c << ") of " << *this);
	return (*get_connection(r, c)).dValue;
}


template<typename T>
bool SparseMatrix<T>::set_as_copy_of(const SparseMatrix<T> &B, double scale)
{

	size_t total=0;
	for(size_t i=0; i<B.rows; i++)
		total += B.num_connections(i);

	if(!is_finalized() || total != consmemsize)
	{
		destroy();
		rows = B.rows;
		cols = B.cols;

		pRowStart = new connection*[rows+1];
		UG_ASSERT(pRowStart != NULL, "out of memory, no more space for " << sizeof(connection*)*(rows+1));

		iMaxNrOfConnections = new size_t[rows];
		UG_ASSERT(iMaxNrOfConnections != NULL, "out of memory, no more space for " << sizeof(size_t)*rows);

		consmem = new connection[total];
		UG_ASSERT(consmem != NULL, "out of memory, no more space for " << sizeof(connection)*total);
	}

	rows = B.rows;
	cols = B.cols;

	connection *p = consmem;

	if(scale == 1.0)
	{
		for(size_t i=0; i<rows; i++)
		{
			size_t nr=B.num_connections(i);
			pRowStart[i] = p;
			iMaxNrOfConnections[i] = nr;

			for(size_t k=0; k < nr; k++, p++)
				*p = B.pRowStart[i][k];

		}
	}
	else
	{
		for(size_t i=0; i<rows; i++)
		{
			size_t nr=B.num_connections(i);
			pRowStart[i] = p;
			iMaxNrOfConnections[i] = nr;

			for(size_t k=0; k < nr; k++, p++)
			{
				*p = B.pRowStart[i][k];
				p->dValue *= scale;
			}

		}
	}
	assert(p == consmem+total);

	pRowStart[rows] = p;
	pRowEnd = pRowStart+1;

	iFragmentedMem = 0;
	iTotalNrOfConnections = total;
	consmemsize = total;


	return true;
}

template<typename T>
bool SparseMatrix<T>::scale(double d)
{
	for(size_t i=0; i < rows; i++)
		for(row_iterator it_k = begin_row(i); it_k != end_row(i); ++it_k)
			(*it_k).dValue *= d;

	return true;
}


} // namespace ug

#endif

