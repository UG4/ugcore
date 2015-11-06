/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__CPU_ALGEBRA__GPUSparseMatrix_IMPL__
#define __H__UG__CPU_ALGEBRA__GPUSparseMatrix_IMPL__

// creation etc

#include <fstream>
#include <cstring>

#include "lib_algebra/common/operations_vec.h"
#include "common/profiler/profiler.h"
#include "gpusparsematrix.h"
#include <vector>
#include <algorithm>
#include "cuda/cuda_manager.h"


namespace ug{

template<typename T>
GPUSparseMatrix<T>::GPUSparseMatrix()
{
	initGPU();
	PROFILE_GPUMATRIX(GPUSparseMatrix_constructor);
	bNeedsValues = true;
	iIterators=0;
	nnz = 0;
	m_numCols = 0;
	maxValues = 0;
	cols.resize(32);
	if(bNeedsValues) values.resize(32);
}

template<typename T>
bool GPUSparseMatrix<T>::resize_and_clear(size_t newRows, size_t newCols)
{
	PROFILE_GPUMATRIX(GPUSparseMatrix_resize_and_clear);
	UG_LOG(this << "GPUSparseMatrix::resize_and_clear(" << newRows << ", " << newCols << ")");
	rowStart.clear(); rowStart.resize(newRows+1, -1);
	rowMax.clear(); rowMax.resize(newRows);
	rowEnd.clear(); rowEnd.resize(newRows, -1);
	m_numCols = newCols;
	nnz = 0;

	cols.clear(); cols.resize(newRows);
	values.clear();
	if(bNeedsValues) values.resize(newRows);
	maxValues = 0;
	return true;
}

template<typename T>
bool GPUSparseMatrix<T>::resize_and_keep_values(size_t newRows, size_t newCols)
{
	PROFILE_GPUMATRIX(GPUSparseMatrix_resize_and_keep_values);
	UG_LOG(this << "GPUSparseMatrix_resize_and_keep_values(" << newRows << ", " << newCols << ")");
	//UG_LOG("GPUSparseMatrix resize " << newRows << "x" << newCols << "\n");
	if(newRows == 0 && newCols == 0)
		return resize_and_clear(0,0);

	if(newRows != num_rows())
	{
		size_t oldrows = num_rows();
		rowStart.resize(newRows+1, -1);
		rowMax.resize(newRows);
		rowEnd.resize(newRows, -1);
		for(size_t i=oldrows; i<newRows; i++)
		{
			rowStart[i] = -1;
			rowEnd[i] = -1;
		}
	}
	if((int)newCols < m_numCols)
		copyToNewSize(get_nnz_max_cols(newCols), newCols);

	m_numCols = newCols;
	return true;
}


template<typename T>
bool GPUSparseMatrix<T>::set_as_transpose_of(const GPUSparseMatrix<value_type> &B, double scale)
{
	PROFILE_GPUMATRIX(GPUSparseMatrix_set_as_transpose_of);
	resize_and_clear(B.num_cols(), B.num_rows());
	/*rowStart.resize(B.num_cols(), 0);
	rowMax.resize(B.num_cols(), 0);
	rowEnd.resize(B.num_cols(), 0);
	nnz = B.nnz;
	bNeedsValues = B.bNeedsValues;
	if(bNeedsValues) values.resize(nnz);
	cols.resize(nnz);
	maxValues = B.maxValues;
	fragmented = 0;
	m_numCols = B.num_rows();
	size_t r, c;
	for(r=0; r<num_rows(); r++)
		for(const_row_iterator it = B.begin_row(r); it != B.end_row(r); ++it)
			rowMax[it.index()]++;
	rowEnd[0] = rowMax[0];
	rowEnd[0] = 0;
	for(c=1; c<B.num_cols(); c++)
	{
		rowStart[c] = rowMax[c-1];
		rowMax[c] = rowStart[c]+rowMax[c];
		rowEnd[c] = rowStart[c];
	}*/

	for(size_t r=0; r<B.num_rows(); r++)
		for(const_row_iterator it = B.begin_row(r); it != B.end_row(r); ++it)
			operator()(it.index(), r) = scale*it.value();
	// todo: sort rows
	return true;
}

template<typename T>
template<typename vector_t>
inline void GPUSparseMatrix<T>::mat_mult_add_row(size_t row, typename vector_t::value_type &dest, double alpha, const vector_t &v) const
{
	for(const_row_iterator conn = begin_row(row); conn != end_row(row); ++conn)
		MatMultAdd(dest, 1.0, dest, alpha, conn.value(), v[conn.index()]);
}

//! calculate x = alpha*x + beta*A*y (A = this matrix)
template<typename T>
template<typename vector_t>
void GPUSparseMatrix<T>::axpy(double alpha, vector_t &x, double beta, const vector_t &y) const
{
	check_device();
	cusparseDcsrmv(CUDAManager::get_cusparseHandle(), CUSPARSE_OPERATION_NON_TRANSPOSE, num_rows(), num_cols(), get_nnz(),
	        &beta, get_matrix_descr(), get_device_value_ptr(), get_device_rowStart(), get_device_cols(), y.get_dev_ptr(),
	            &alpha, x.get_dev_ptr());
}


// calculate dest = alpha1*v1 + beta1*A*w1 (A = this matrix)
template<typename T>
template<typename vector_t>
bool GPUSparseMatrix<T>::axpy(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const vector_t &w1) const
{
//	UG_LOG("GPUSparseMatrix axpy\n");
	PROFILE_GPUMATRIX(GPUSparseMatrix_axpy);

	if(alpha1 == 0.0 || &dest == &v1)
		axpy(alpha1, dest, beta1, w1);
	else
	{
		axpy(0.0, dest, beta1, w1);
		dest.add(alpha1, v1);
	}
	return true;
}

// calculate dest = alpha1*v1 + beta1*A^T*w1 (A = this matrix)
template<typename T>
template<typename vector_t>
bool GPUSparseMatrix<T>::axpy_transposed(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const vector_t &w1) const
{
	PROFILE_GPUMATRIX(GPUSparseMatrix_axpy_transposed);
	UG_ASSERT(0, "not implemented");
	return true;
}

template<typename T>
bool GPUSparseMatrix<T>::set(double a)
{
	PROFILE_GPUMATRIX(GPUSparseMatrix_set);
	check_fragmentation();
	for(size_t row=0; row<num_rows(); row++)
		for(row_iterator it = begin_row(row); it != end_row(row); ++it)
		{
			if(it.index() == row)
				it.value() = a;
			else
				it.value() = 0.0;
		}

	return true;
}


template<typename T>
inline bool GPUSparseMatrix<T>::is_isolated(size_t i) const
{
	UG_ASSERT(i < num_rows() && i >= 0, *this << ": " << i << " out of bounds.");

	for(const_row_iterator it = begin_row(i); it != end_row(i); ++it)
		if(it.index() != i && it.value() != 0.0)
			return false;
	return true;
}


template<typename T>
void GPUSparseMatrix<T>::set_matrix_row(size_t row, connection *c, size_t nr)
{
	/*PROFILE_GPUMATRIX(GPUSparseMatrix_set_matrix_row);
	bool bSorted=false;
	for(size_t i=0; bSorted && i+1 < nr; i++)
		bSorted = c[i].iIndex < c[i+1].iIndex;
	if(bSorted)
	{
		int start;
		if(rowStart[row] == -1 || rowMax[row] - rowStart[row] < (int)nr)
		{
			assureValuesSize(maxValues+nr);
			start = maxValues;
			rowMax[row] = start+nr;
			rowStart[row] = start;
			maxValues+=nr;
		}
		else
			start = rowStart[row];
		rowEnd[row] = start+nr;
		for(size_t i=0; i<nr; i++)
		{
			cols[start+i] = c[i].iIndex;
			values[start+i] = c[i].dValue;
		}
	}
	else*/
	{
		for(size_t i=0; i<nr; i++)
			operator()(row, c[i].iIndex) = c[i].dValue;
	}
}

template<typename T>
void GPUSparseMatrix<T>::add_matrix_row(size_t row, connection *c, size_t nr)
{
	//PROFILE_GPUMATRIX(GPUSparseMatrix_add_matrix_row);
	for(size_t i=0; i<nr; i++)
		operator()(row, c[i].iIndex) += c[i].dValue;
}


template<typename T>
bool GPUSparseMatrix<T>::set_as_copy_of(const GPUSparseMatrix<T> &B, double scale)
{
	//PROFILE_GPUMATRIX(GPUSparseMatrix_set_as_copy_of);
	resize_and_clear(B.num_rows(), B.num_cols());
	for(size_t i=0; i < B.num_rows(); i++)
		for(const_row_iterator it = B.begin_row(i); it != B.end_row(i); ++it)
			operator()(i, it.index()) = scale*it.value();
	return true;
}



template<typename T>
bool GPUSparseMatrix<T>::scale(double d)
{
	//PROFILE_GPUMATRIX(GPUSparseMatrix_scale);
	for(size_t i=0; i < num_rows(); i++)
		for(row_iterator it = begin_row(i); it != end_row(i); ++it)
			it.value() *= d;
	return true;
}





//======================================================================================================
// submatrix set/get

template<typename T>
template<typename M>
void GPUSparseMatrix<T>::add(const M &mat)
{
	for(size_t i=0; i < mat.num_rows(); i++)
	{
		int r = mat.row_index(i);
		for(size_t j=0; j < mat.num_cols(); j++)
		{
			int c = mat.col_index(j);
			(*this)(r,c) += mat(i,j);
		}
	}
}


template<typename T>
template<typename M>
void GPUSparseMatrix<T>::set(const M &mat)
{
	for(size_t i=0; i < mat.num_rows(); i++)
	{
		int r = mat.row_index(i);
		for(size_t j=0; j < mat.num_cols(); j++)
		{
			int c = mat.col_index(j);
			(*this)(r,c) = mat(i,j);
		}
	}
}

template<typename T>
template<typename M>
void GPUSparseMatrix<T>::get(M &mat) const
{
	for(size_t i=0; i < mat.num_rows(); i++)
	{
		int r = mat.row_index(i);
		for(size_t j=0; j < mat.num_cols(); j++)
		{
			int c = mat.col_index(j);
			mat(i,j) = (*this)(r,c);
		}
	}
}


template<typename T>
int GPUSparseMatrix<T>::get_index_internal(size_t row, int col) const
{
//	PROFILE_GPUMATRIX(SP_get_index_internal);
	assert(rowStart[row] != -1);
	int l = rowStart[row];
	int r = rowEnd[row];
	int mid=0;
	while(l < r)
	{
		mid = (l+r)/2;
		if(cols[mid] < col)
			l = mid+1;
		else if(cols[mid] > col)
			r = mid-1;
		else
			return mid;
	}
	mid = (l+r)/2;
	int ret;
	if(mid < rowStart[row])
		ret = rowStart[row];
	else if(mid == rowEnd[row] || col <= cols[mid])
		ret = mid;
	else ret = mid+1;
	UG_ASSERT(ret <= rowEnd[row] && ret >= rowStart[row], "row iterator row " <<  row << " pos " << ret << " out of bounds [" << rowStart[row] << ", " << rowEnd[row] << "]");
	return ret;
}



template<typename T>
int GPUSparseMatrix<T>::get_index_const(int r, int c) const
{
	if(rowStart[r] == -1 || rowStart[r] == rowEnd[r]) return -1;
	int index=get_index_internal(r, c);
	if(index >= rowStart[r] && index < rowEnd[r] && cols[index] == c)
		return index;
	else
		return -1;
}


template<typename T>
int GPUSparseMatrix<T>::get_index(int r, int c)
{
//	UG_LOG("get_index " << r << ", " << c << "\n");
//	UG_LOG(rowStart[r] << " - " << rowMax[r] << " - " << rowEnd[r] << " - " << cols.size() << " - "  << maxValues << "\n");
	if(rowStart[r] == -1 || rowStart[r] == rowEnd[r])
	{
//		UG_LOG("new row\n");
		// row did not start, start new row at the end of cols array
		assureValuesSize(maxValues+1);
		rowStart[r] = maxValues;
		rowEnd[r] = maxValues+1;
		rowMax[r] = maxValues+1;
		if(bNeedsValues) values[maxValues] = 0.0;
		cols[maxValues] = c;
		maxValues++;
		nnz++;
		return maxValues-1;
	}

	/*    for(int i=rowStart[r]; i<rowEnd[r]; i++)
	 if(cols[i] == c)
	 return i;*/

	// get the index where (r,c) _should_ be
	int index=get_index_internal(r, c);

	if(index < rowEnd[r]
			&& cols[index] == c)
	{
//		UG_LOG("found\n");
		return index; // we found it
	}

	// we did not find it, so we have to add it

#ifndef NDEBUG
	assert(index == rowEnd[r] || cols[index] > c);
	assert(index == rowStart[r] || cols[index-1] < c);
	for(int i=rowStart[r]+1; i<rowEnd[r]; i++)
		assert(cols[i] > cols[i-1]);
#endif
	if(rowEnd[r] == rowMax[r] && rowEnd[r] == maxValues
			&& maxValues < (int)cols.size())
	{
//		UG_LOG("at end\n");
		// this row is stored at the end, so we can easily make more room
		rowMax[r]++;
		maxValues++;
	}

	if(rowEnd[r] == rowMax[r])
	{
//		UG_LOG("renew\n");
		// row is full, we need to copy it to another place
		int newSize = (rowEnd[r]-rowStart[r])*2;
		if(maxValues+newSize > (int)cols.size())
		{
//			UG_LOG("ass\n");
			assureValuesSize(maxValues+newSize);
			index=get_index_internal(r, c);
		}
		// copy it to the end and insert index
		fragmented += rowEnd[r]-rowStart[r];
		index = index-rowStart[r]+maxValues;
		int j=rowEnd[r]-rowStart[r]+maxValues;
		if(rowEnd[r] != 0)
			for(int i=rowEnd[r]-1; i>=rowStart[r]; i--, j--)
			{
				if(j==index) j--;
				if(bNeedsValues) values[j] = values[i];
				cols[j] = cols[i];
				if(i==rowStart[r]) break;
			}
		rowEnd[r] = maxValues+rowEnd[r]-rowStart[r]+1;
		rowStart[r] = maxValues;
		rowMax[r] = maxValues+newSize;
		maxValues += newSize;
	}
	else
	{
//		UG_LOG("enlength\n");
		// move part > index so we can insert the index
		if(rowEnd[r] != 0)
			for(int i=rowEnd[r]-1; i>=index; i--)
			{
				if(bNeedsValues) values[i+1] = values[i];
				cols[i+1] = cols[i];
				if(i==index) break;
			}
		rowEnd[r]++;
	}
	if(bNeedsValues) values[index] = 0.0;
	cols[index] = c;

	nnz++;
#ifndef NDEBUG
	assert(index >= rowStart[r] && index < rowEnd[r]);
	for(int i=rowStart[r]+1; i<rowEnd[r]; i++)
		assert(cols[i] > cols[i-1]);
#endif
	return index;

}

template<typename T>
void GPUSparseMatrix<T>::copyToNewSize(size_t newSize, size_t maxCol)
{
	PROFILE_GPUMATRIX(GPUSparseMatrix_copyToNewSize);
	/*UG_LOG("copyToNewSize: from " << values.size()  << " to " << newSize << "\n");
	UG_LOG("sizes are " << cols.size() << " and " << values.size() << ", ");
	UG_LOG(reset_floats << "capacities are " << cols.capacity() << " and " << values.capacity() << ", NNZ = " << nnz << ", fragmentation = " <<
			(1-((double)nnz)/cols.size())*100.0 << "%\n");
	if(newSize == nnz) { UG_LOG("Defragmenting to NNZ."); }*/
	if( (iIterators > 0)
		|| (newSize > values.size() && (100.0*nnz)/newSize < 20 && newSize <= cols.capacity()) )
	{
		UG_ASSERT(newSize > values.size(), "no nnz-defragmenting while using iterators.");
		//UG_LOG("copyToNew not defragmenting because of iterators or low fragmentation.\n");}
		cols.resize(newSize);
		cols.resize(cols.capacity());
		if(bNeedsValues) { values.resize(newSize); values.resize(cols.size()); }
		return;
	}

	std::vector<value_type> v(newSize);
	std::vector<int> c(newSize);
	size_t j=0;
	for(size_t r=0; r<num_rows(); r++)
	{
		if(rowStart[r] == -1)
			rowStart[r] = rowEnd[r] = rowMax[r] = j;
		else
		{
			size_t start=j;
			for(int k=rowStart[r]; k<rowEnd[r]; k++)
			{
				if(cols[k] < (int)maxCol)
				{
					if(bNeedsValues) v[j] = values[k];
					c[j] = cols[k];
					j++;
				}
			}
			rowStart[r] = start;
			rowEnd[r] = rowMax[r] = j;
		}
	}
	rowStart[num_rows()] = rowEnd[num_rows()-1];
	fragmented = 0;
	maxValues = j;
	if(bNeedsValues) std::swap(values, v);
	std::swap(cols, c);
}

template<typename T>
void GPUSparseMatrix<T>::check_fragmentation() const
{
	if((double)nnz/(double)maxValues < 0.9)
		(const_cast<this_type*>(this))->defragment();
}

template<typename T>
void GPUSparseMatrix<T>::assureValuesSize(size_t s)
{
    if(s <= cols.size()) return;
    size_t newSize = nnz*2;
    if(newSize < s) newSize = s;
    copyToNewSize(newSize);

}

template<typename T>
int GPUSparseMatrix<T>::get_nnz_max_cols(size_t maxCols)
{
	int j=0;
	for(size_t r=0; r<num_rows(); r++)
	{
		if(rowStart[r] == -1) continue;
		for(int k=rowStart[r]; k<rowEnd[r]; k++)
			if(cols[k] < (int)maxCols)
				j++;
	}
	return j;
}

} // namespace ug

#endif

