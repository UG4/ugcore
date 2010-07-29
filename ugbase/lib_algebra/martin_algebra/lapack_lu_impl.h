/*
 *  lapack_lu_impl.h
 *
 *  Created by Martin Rupp on 26.11.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 */

#ifndef __H__UG__MARTIN_ALGEBRA__LAPACK_LU_IMPL__
#define __H__UG__MARTIN_ALGEBRA__LAPACK_LU_IMPL__



//using namespace std;
namespace ug{

extern "C"
{
    void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
    void dgetrs_(char *trans, int *n, int *nrhs, const double *a, int *lda,
            const int *ipiv, double *b, int *ldb, int *info);
}


template<typename vec_type>
void LapackLU::apply(const vec_type &b, vec_type &x)
{
#ifndef NDEBUG
    const size_t nrOfUnknowns = block_vector_traits<typename vec_type::entry_type>::nrOfUnknowns;
	UG_ASSERT(size == b.size() * nrOfUnknowns && size == x.size() * nrOfUnknowns, " wrong size! has to be " << size << ", but is " << b << " and " << x);
#endif

	x = b;
	// TODO: this only works for fixed array entries.

	// solve system
	char trans ='N';
	int dim = size;
	int nrhs = 1;
	int info;

	dgetrs_(&trans, &dim, &nrhs, densemat, &dim, interchange, (double*)&x[0], &dim, &info);
	UG_ASSERT(info == 0, "info is " << info);
}


template<typename matrix_type>
void LapackLU::init(const matrix_type &A)
{
  const size_t nrOfUnknowns = block_matrix_traits<typename matrix_type::entry_type>::nrOfUnknowns;
	size = A.num_rows() * nrOfUnknowns;

	if(densemat) delete[] densemat;
	if(interchange) delete[] interchange;

	densemat = new double[size*size];
	interchange = new int[size];

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

	int info = 0;
	int dim = size;
	dgetrf_(&dim, &dim, densemat, &dim, interchange, &info);
	UG_ASSERT(info == 0, "info is " << info << ( info > 0 ? ": matrix singular in U(i,i)" : ": i-th argument had had illegal value"));
}

} // namespace ug

#endif /* __H__UG__MARTIN_ALGEBRA__LAPACK_LU_IMPL__ */
