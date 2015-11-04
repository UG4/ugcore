/**
 * \file sparsematrix_print.h
 *
 * \author Martin Rupp
 *
 * \date 29.10.2012
 *
 * Goethe-Center for Scientific Computing 2012
 */

#ifndef __H__UG__CPU_ALGEBRA__SPARSEMATRIX_PRINT__
#define  __H__UG__CPU_ALGEBRA__SPARSEMATRIX_PRINT__

#include "sparsematrix.h"
#include "common/common.h"

namespace ug {

/// \addtogroup crs_algebra
/// \{

//!
//! print to console whole SparseMatrix
template<typename T>
void SparseMatrix<T>::print(const char * const text) const
{
	UG_LOG("================= SparseMatrix " << num_rows() << "x" << num_cols() << " =================\n");
	for(size_t i=0; i < num_rows(); i++)
		printrow(i);
}


//!
//! print the row row to the console
template<typename T>
void SparseMatrix<T>::printrow(size_t row) const
{
	UG_LOG("row " << row << ": ");
	for(const_row_iterator it=begin_row(row); it != end_row(row); ++it)
	{
		if(it.value() == 0.0) continue;
		UG_LOG(" ");
		UG_LOG("(" << it.index() << " -> " << it.value() << ")");
	}

	UG_LOG("\n");
}

template<typename T>
void SparseMatrix<T>::printtype() const
{
	std::cout << *this;
}

// end group crs_algebra
/// \}

}
#endif // __H__UG__CPU_ALGEBRA__SPARSEMATRIX_PRINT__
