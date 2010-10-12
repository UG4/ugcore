/**
 * \file lapack_lu_impl.h
 *
 * \author Martin Rupp
 *
 * \date 26.11.2009
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#ifndef __H__UG__MARTIN_ALGEBRA__LAPACK_LU_IMPL__
#define __H__UG__MARTIN_ALGEBRA__LAPACK_LU_IMPL__

#include "lapack/lapack.h"

//using namespace std;
namespace ug{



template<typename vec_type>
bool LapackLU::apply(const vec_type &b, vec_type &x)
{
#ifndef NDEBUG
    const size_t nrOfUnknowns = block_vector_traits<typename vec_type::entry_type>::nrOfUnknowns;
	UG_ASSERT(size == b.size() * nrOfUnknowns && size == x.size() * nrOfUnknowns, " wrong size! has to be " << size << ", but is " << b << " and " << x);
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

	UNUSED_VARIABLE(info);

//	we're done
	return true;
}


template<typename matrix_type>
bool LapackLU::init(const matrix_type &A)
{
	// TODO: Use nrOfRows and nrOfCols here. I have changed this due to compile errors. Andreas Vogel
	// before: const size_t nrOfUnknowns = block_matrix_traits<typename matrix_type::entry_type>::nrOfUnknowns;
	const size_t nrOfUnknowns = block_matrix_traits<typename matrix_type::entry_type>::nrOfRows;
	size = A.num_rows() * nrOfUnknowns;

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
			int rr = r*nrOfUnknowns;
			int cc = (*it).iIndex*nrOfUnknowns;
			for(size_t r2=0; r2<nrOfUnknowns; r2++)
					for(size_t c2=0; c2<nrOfUnknowns; c2++)
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

#endif /* __H__UG__MARTIN_ALGEBRA__LAPACK_LU_IMPL__ */
