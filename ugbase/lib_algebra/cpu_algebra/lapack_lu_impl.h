/**
 * \file lapack_lu_impl.h
 *
 * \author Martin Rupp
 *
 * \date 26.11.2009
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#ifndef __H__UG__CPU_ALGEBRA__LAPACK_LU_IMPL__
#define __H__UG__CPU_ALGEBRA__LAPACK_LU_IMPL__

#include "../small_algebra/lapack/lapack.h"

//using namespace std;
namespace ug{



template<typename vec_type>
bool LapackLU::apply(const vec_type &b, vec_type &x)
{
#ifndef NDEBUG
	const size_t static_size = block_vector_traits<typename vec_type::value_type>::static_size;
	UG_ASSERT(size == b.size() * static_size && size == x.size() * static_size,
			" wrong size! has to be " << size << ", but is " << b << " and " << x);
#endif

	x = b;
	// TODO: this only works for fixed array entries.
	int info = getrs(ModeNoTrans, size, 1, densemat, size, interchange, (double*)&x[0], size);
	UG_ASSERT(info == 0, "info is " << info);
	if(info != 0)
	{
		UG_LOG("ERROR in LapackLU::apply: Cannot apply lapack routine.\n");
		return false;
	}

	//UNUSED_VARIABLE(info);

//	we're done
	return true;
}


template<typename matrix_type>
bool LapackLU::init(const matrix_type &A)
{
	UG_ASSERT(block_matrix_traits<typename matrix_type::value_type>::is_static,
			"non-static blockmatrices not fully implemented yet\n");
	const size_t nrOfRows = block_matrix_traits<typename matrix_type::value_type>::static_num_rows;
	UG_ASSERT(nrOfRows == block_matrix_traits<typename matrix_type::value_type>::static_num_cols, "only square matrices supported");
	size = A.num_rows() * nrOfRows;

	/*else
	{
		UG_ASSERT(0,
		for(size_t r=0; r<A.num_rows(); r++)
		{
			typename matrix_type::cRowIterator it(A, r);
			if(it.isEnd()) continue;
			size += GetRows((*it).dValue);
		}
	}*/

	if(densemat) delete[] densemat;
	if(interchange) delete[] interchange;

	densemat = new double[size*size];
	interchange = new int[size];
	if(densemat == NULL || interchange == NULL)
	{
		UG_LOG("ERROR in LapackLU::init: Not enough memory to allocate lapack matrices.\n");
		return false;
	}

	memset(densemat, 0, sizeof(double)*size*size);

	for(size_t r=0; r<A.num_rows(); r++)
		for(typename matrix_type::cRowIterator it(A, r); !it.isEnd(); ++it)
		{
			int rr = r*nrOfRows;
			int cc = (*it).iIndex*nrOfRows;
			for(size_t r2=0; r2<nrOfRows; r2++)
					for(size_t c2=0; c2<nrOfRows; c2++)
					  densemat[ (rr + r2) + (cc+c2)*size ] = BlockRef((*it).dValue, r2, c2);
		}

	int info = getrf(size, size, densemat, size, interchange);
	UG_ASSERT(info == 0, "info is " << info << ( info > 0 ? ": matrix singular in U(i,i)" : ": i-th argument had an illegal value"));
	// todo: put this error text in a lapack_errors.cpp file
	if(info != 0)
	{
		UG_LOG("ERROR in LapackLU::init: Cannot init lapack routine.\n");
		return false;
	}

	UNUSED_VARIABLE(info);

//	we're done
	return true;
}

} // namespace ug

#endif /* __H__UG__CPU_ALGEBRA__LAPACK_LU_IMPL__ */
