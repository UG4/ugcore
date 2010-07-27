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


//template<typename vec_type>
void LapackLU::apply(const Vector<double> &b, Vector<double> &x)
{
	//cout << "LapackLU::apply" << endl;
	// TODO: Variing nr of unknowns
	//int nrOfUnknowns = block_vector_traits<vec_type>::nrOfUnknowns;
#ifndef NDEBUG
	const int nrOfUnknowns = 1;
	UG_ASSERT(size == b.size() * nrOfUnknowns && size == x.size() * nrOfUnknowns, " wrong size! has to be " << size << ", but is " << b << " and " << x);
#endif
	x = b;
	//for(int i=0; i < b.size(); i++)
		//for(int j=0; j<nrOfUnknowns; j++)
			//vec[i*nrOfUnknowns+j] = getAt(b[i], j);

	// solve system
	char trans ='N';
	int dim = size;
	int nrhs = 1;
	int info;

	dgetrs_(&trans, &dim, &nrhs, densemat, &dim, interchange, &x[0], &dim, &info);
	UG_ASSERT(info == 0, "info is " << info);
	//for(int i=0; i<x.size(); i++)
		//for(int j=0; j<nrOfUnknowns; j++)
			//setAt(x[i], j, vec[i*nrOfUnknowns+j]);

}


//template<typename entry_type>
void LapackLU::init(const SparseMatrix<double> &A)
{
	//cout << "LapackLU::init, A = " << A.num_rows() << " x " << A.num_cols() << endl;
	const int nrOfUnknowns = 1 ; //block_matrix_traits<entry_type>::nrOfUnknowns;
	size = A.num_rows() * nrOfUnknowns;

	if(densemat) delete[] densemat;
	if(interchange) delete[] interchange;

	densemat = new double[size*size];
	interchange = new int[size];

	memset(densemat, 0, sizeof(double)*size*size);

	/*for(int r=0; r<A.num_rows(); r++)
		for(typename SparseMatrix<entry_type>::cRowIterator it(A, r); !it.isEnd(); ++it)
		{
			int rr = r*nrOfUnknowns;
			int cc = (*it).iIndex*nrOfUnknowns;
			for(int r2=0; r2<nrOfUnknowns; r2++)
					for(int c2=0; c2<nrOfUnknowns; c2++)
						densemat[(rr+r2) + (cc+c2)*size] = getAt((*it).dValue, r2, c2);
		}*/

	for(size_t r=0; r<A.num_rows(); r++)
		for(SparseMatrix<double>::cRowIterator it = A.beginRow(r); !it.isEnd(); ++it)
		{
			int rr = r*nrOfUnknowns;
			int cc = (*it).iIndex*nrOfUnknowns;
			densemat[rr + cc*size] = (*it).dValue;
		}

	int info = 0;
	int dim = size;
	dgetrf_(&dim, &dim, densemat, &dim, interchange, &info);
	UG_ASSERT(info == 0, "info is " << info << ( info > 0 ? ": matrix singular in U(i,i)" : ": i-th argument had had illegal value"));
}

} // namespace ug

#endif /* __H__UG__MARTIN_ALGEBRA__LAPACK_LU_IMPL__ */
