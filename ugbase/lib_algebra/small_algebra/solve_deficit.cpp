#include "solve_deficit.h"

#include "common/log.h"
#include "common/debug_print.h"
#include "lib_algebra/common/operations_vec.h"

typedef int lapack_int;
typedef float lapack_float;
typedef double lapack_double;
typedef int lapack_ftnlen;


inline double dabs(double a) { return a > 0 ? a : -a; }
/*
extern "C"
{

	// factor system *GETRF
void dsytrs_(char *uplo, int *n, int *nrhs, double *a, int *lda,
             int *ipivot, double *b, int *ldb, int *info);
}


inline lapack_int sytrs(char uplo, int n, int nrhs, double *a, int  lda,
        int *ipivot, double *b, int ldb)
{
    lapack_int info;
    dsytrs_(&uplo, &n, &nrhs, a, &lda, ipivot, b, &ldb, &info);
    return info;
}

bool SolveDeficit2(DenseMatrix< VariableArray2<double> > &A,
		DenseVector<VariableArray1<double> > &x, DenseVector<VariableArray1<double> > &rhs)
{
	DenseMatrix< VariableArray2<double> > A2=A;
	DenseVector<VariableArray1<double> > rhs2=rhs;

	std::vector<int> interchange(A.num_rows());
	int info = sytrs('U', A.num_rows(), 1, &A(0,0), A.num_rows(), &interchange[0], &rhs[0], A.num_rows());
	x=rhs;
	DenseVector<VariableArray1<double> > f;
		f = A2*x - rhs2;
		if(VecNormSquared(f) > 1e-12)
		{
			UG_LOG("solving was wrong:");
			UG_LOGN(CPPString(A2, "Aold"));
			rhs2.maple_print("rhs");
			UG_LOGN(CPPString(A, "Adecomp"));
			rhs.maple_print("rhsDecomp");
			x.maple_print("x");
			f.maple_print("f");

		}

}*/

namespace ug{

/*
|   free variable
xxxxxxxx
0xxxxxxx
000xxxxx
0000xxxx
00000xxx	<- iNonNullRows = 5
00000000 	\
00000000	|	lin. dep.
00000000	/

	 iNonNullRows = number of bounded variable (that is: number of linear independent rows)
	 returns true if solveable.

 */
bool Decomp(DenseMatrix< VariableArray2<double> > &A, DenseVector<VariableArray1<double> > &rhs,
		size_t &iNonNullRows, std::vector<size_t> &xinterchange)
{
	// LU Decomposition, IKJ Variant
	size_t rows = A.num_rows();
	size_t cols = A.num_cols();

	xinterchange.resize(cols);
	for(size_t k=0; k<cols; k++)
		xinterchange[k] = k;

//	UG_LOG("DECOMP " << rows << " x " << cols << "\n");

	std::vector<size_t> colInterchange;
	size_t k;
	for(k=0; k<cols; k++)
	{
//		UG_LOG("-----------------" << k << " < " << cols << " ------\n");
		size_t biggestR, biggestC;
		double dBiggest = -1;

		for(size_t r=k; r<rows; r++)
			for(size_t c=k; c<cols; c++)
			{
				double n = dabs(A(r,c));
				if(dBiggest < n)
				{
					dBiggest=n;
					biggestR = r;
					biggestC = c;
				}
			}


//		UG_LOG(CPPString(A, "BeforeSwap"));
		if(dBiggest < 1e-14)
		{
//			UG_LOG("k = " << k << " abort");
			break;
		}
//		UG_LOG("k = " << k << " biggest = " << biggestR << ", " << biggestC << "\n");
		if(biggestR != k)
		{
			for(size_t j=0; j<cols; j++)
				std::swap(A(k, j), A(biggestR, j));
			std::swap(rhs[k], rhs[biggestR]);
		}
		if(biggestC != k)
		{
			for(size_t j=0; j<rows; j++)
				std::swap(A(j, k), A(j, biggestC));
			std::swap(xinterchange[k], xinterchange[biggestC]);
		}
//		UG_LOG(CPPString(A, "AfterSwap"));

		for(size_t i=k+1; i<rows; i++)
		{
			double a = A(i,k)/A(k,k);
			for(size_t j=k+1; j<cols; j++)
				A(i,j) = A(i,j) - a*A(k,j);
			rhs[i] -= a*rhs[k];
			A(i,k) = 0;

		}


//		UG_LOG(CPPString(A, "MAT"));
//		PrintVector(rhs, "rhs");
	}
	iNonNullRows = k;
	// every row below iNonNullRows is zero.
	// equation system is solveable if every rhs[k] = 0 for k>=iNonNullRows
	for(; k<rows; k++)
		if(dabs(rhs[k]) > 1e-7 )
			return false;
		else
			rhs[k] = 0.0;
//	UG_LOGN("iNonNullRows = " << iNonNullRows << " abort");
	return true;
}


bool SolveDeficit(DenseMatrix< VariableArray2<double> > &A,
		DenseVector<VariableArray1<double> > &x, DenseVector<VariableArray1<double> > &rhs)
{
	DenseMatrix< VariableArray2<double> > A2=A;
	DenseVector<VariableArray1<double> > rhs2=rhs;

	UG_ASSERT(A.num_rows() == rhs.size(), "");
	UG_ASSERT(A.num_cols() == x.size(), "");

	size_t iNonNullRows;
	x.resize(A.num_cols());
	for(size_t i=0; i<x.size(); i++)
		x[i] = 0.0;
	std::vector<size_t> interchange;
	if(Decomp(A, rhs, iNonNullRows, interchange) == false) return false;

//	A.maple_print("Adecomp");
//	rhs.maple_print("rhs decomp");

	for(int i=iNonNullRows-1; i>=0; i--)
	{
		double d=A(i,i);
		double s=0;
		for(int k=i+1; k<A.num_cols(); k++)
			s += A(i,k)*x[interchange[k]];
		x[interchange[i]] = (rhs[i] - s)/d;
	}
	DenseVector<VariableArray1<double> > f;
	f = A2*x - rhs2;
	if(VecNormSquared(f) > 1e-12)
	{
		UG_LOGN("iNonNullRows = " << iNonNullRows);
		UG_LOG("solving was wrong:");
		UG_LOGN(CPPString(A2, "Aold"));
		rhs2.maple_print("rhs");
		UG_LOGN(CPPString(A, "Adecomp"));
		rhs.maple_print("rhsDecomp");
		x.maple_print("x");
		f.maple_print("f");

	}

	return true;
}
}
