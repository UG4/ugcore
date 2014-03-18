#ifndef __H__UG__CPU_ALGEBRA__LU_DECOMP_H__
#define __H__UG__CPU_ALGEBRA__LU_DECOMP_H__

#include "../small_matrix/densematrix.h"
#include "../small_matrix/densevector.h"
#include "../../common/operations.h"

#include <vector>

#ifndef LAPACK_AVAILABLE

namespace ug
{

template<typename matrix_t>
bool LUDecomp(DenseMatrix<matrix_t> &A, size_t *pInterchange)
{
	// LU Decomposition, IKJ Variant
	UG_ASSERT(A.num_rows() == A.num_cols(), "LU decomposition only for square matrices");
	if(A.num_rows() != A.num_cols()) return false;

	size_t n = A.num_rows();

	if(pInterchange)
	{
		pInterchange[0] = 0;
		for(size_t k=0; k<n; k++)
		{
			size_t biggest = k;
			for(size_t j=k+1; j<n; j++)
				if(dabs(A(biggest, k)) < dabs(A(j,k)))
					biggest = j; // costly.
			if(biggest != k)
				for(size_t j=0; j<n; j++)
					std::swap(A(k, j), A(biggest, j));

			pInterchange[k] = biggest;
			if(dabs(A(k,k)) < 1e-10)
				return false;

			for(size_t i=k+1; i<n; i++)
			{
				A(i,k) = A(i,k)/A(k,k);
				for(size_t j=k+1; j<n; j++)
					A(i,j) = A(i,j) - A(i,k)*A(k,j);
			}
		}
	}
	else
	{
		for(size_t k=0; k<n; k++)
		{
			if(dabs(A(k,k)) < 1e-10)
				return false;

			for(size_t i=k+1; i<n; i++)
			{
				A(i,k) = A(i,k)/A(k,k);
				for(size_t j=k+1; j<n; j++)
					A(i,j) = A(i,j) - A(i,k)*A(k,j);
			}
		}
	}
	return true;
}

template<typename matrix_t>
bool LUDecompIKJ(DenseMatrix<matrix_t> &A, size_t *pInterchange)
{
	// LU Decomposition, IKJ Variant
	UG_ASSERT(A.num_rows() == A.num_cols(), "LU decomposition only for square matrices");
	if(A.num_rows() != A.num_cols()) return false;

	size_t n = A.num_rows();

	if(pInterchange)
	{
		pInterchange[0] = 0;
		for(size_t i=0; i < n; i++)
		{
			UG_LOG("i=" << i << ": \n")
			size_t biggest = i;
			UG_LOG("A(i,i) = " << A(i,i) << "\n");
			for(size_t j=i+1; j<n; j++)
			{
				UG_LOG("A(" << j << ", " << i << ")=" << A(j,i) << "\n");
				if(dabs(A(biggest, i)) < dabs(A(j,i))) biggest = j; // costly.
			}
			UG_LOG(biggest << " is biggest.");

			if(biggest != i)
				for(size_t j=0; j<n; j++)
					std::swap(A(i, j), A(biggest, j));

			pInterchange[i] = biggest;
			if(dabs(A(i,i)) < 1e-10)
				return false;

			// eliminate all entries A(i, k) with k<i with rows A(k, .) and k<i
			for(size_t k=0; k<i; k++)
			{
				// add row k to row i by A(i, .) -= A(k,.)  A(i,k) / A(k,k)
				// so that A(i,k) is zero (elements A(i, <k) are already "zero")
				// safe A(i,k)/A(k,k) in A(i,k)
				double a_ik = (A(i,k) /= A(k,k));

				for(size_t j=k+1; j<n; j++)
					A(i,j) -= A(k,j) * a_ik;
			}
		}
		// P A = L R
	}
	else
	{
		// for all rows
		for(size_t i=0; i < n; i++)
		{
			if(dabs(A(i,i)) < 1e-10)
				return false;
			// eliminate all entries A(i, k) with k<i with rows A(k, .) and k<i
			for(size_t k=0; k<i; k++)
			{
				// add row k to row i by A(i, .) -= A(k,.)  A(i,k) / A(k,k)
				// so that A(i,k) is zero (elements A(i, <k) are already "zero")
				// safe A(i,k)/A(k,k) in A(i,k)
				double a_ik = (A(i,k) /= A(k,k));

				for(size_t j=k+1; j<n; j++)
					A(i,j) -= A(k,j) * a_ik;
			}
		}
	}

	return true;
}
template<typename matrix_t, typename vector_t>
bool SolveLU(const DenseMatrix<matrix_t> &A, DenseVector<vector_t> &x, const size_t *pInterchange)
{
	size_t n=A.num_rows();

	if(pInterchange)
		for(size_t i=0; i<n; i++)
			if(i < pInterchange[i])
				std::swap(x[i], x[pInterchange[i]]);

	// LU x = b, -> U x = L^{-1} b
	// solve lower left
	double s;
	for(size_t i=0; i<n; i++)
	{
		s = x[i];
		for(size_t k=0; k<i; k++)
			s -= A(i, k)*x[k];
		x[i] = s;
	}

	// -> x = U^{-1} (L^{-1} b)
	// solve upper right
	for(size_t i=n-1; ; i--)
	{
		s = x[i];
		for(size_t k=i+1; k<n; k++)
			s -= A(i, k)*x[k];
		x[i] = s/A(i,i);
		if(i==0) break;
	}

	return true;
}


///////////////////////////////////////////////////////////////////////////////////////
/**
 * smallInverse<size_t n>
 * A class to hold a inverse of a smallMatrix<n>
 */
template<typename TStorage>
class DenseMatrixInverse
{
private: // storage
	DenseMatrix<TStorage> densemat;
	std::vector<size_t> interchange;

public:
	inline size_t num_cols() const
	{
		return densemat.num_cols();
	}

	inline size_t num_rows() const
	{
		return densemat.num_rows();
	}

	inline void resize(size_t k, size_t k2)
	{
		UG_COND_THROW(k!=k2, "only square matrices supported");
		resize(k);
	}
	inline void resize(size_t k)
	{
		densemat.resize(k,k);
		densemat = 0.0;
	}

	double &operator()(int r, int c)
	{
		return densemat(r,c);
	}
	const double &operator()(int r, int c) const
	{
		return densemat(r,c);
	}
public:
	//! initializes this object as inverse of mat
	bool set_as_inverse_of(const DenseMatrix<TStorage> &mat)
	{
		UG_ASSERT(mat.num_rows() == mat.num_cols(), "only for square matrices");
		densemat = mat;
		return invert();
	}

	bool invert()
	{
		interchange.resize(densemat.num_rows());

		if(!interchange.empty()){
			bool bLUDecomp = LUDecomp(densemat, &interchange[0]);

			if(bLUDecomp!=true)
			{
				UG_LOG("ERROR in 'DenseMatrixInverse::invert': Matrix is singular, "
						"cannot calculate Inverse.\n");
			}

			return bLUDecomp;
		}
		else
			return true;
	}

	template<typename vector_t>
	void apply(DenseVector<vector_t> &vec) const
	{
		if(!interchange.empty())
			SolveLU(densemat, vec, &interchange[0]);
	}

	// todo: implem

	// todo: implement operator *=

	template<typename T> friend std::ostream &operator << (std::ostream &out, const DenseMatrixInverse<T> &mat);
};


template<typename T>
std::ostream &operator << (std::ostream &out, const DenseMatrixInverse<T> &mat)
{
	out << "[ ";
	for(size_t r=0; r<mat.num_rows(); ++r)
	{
		for(size_t c=0; c<mat.num_cols(); ++c)
			out << mat.densemat(r, c) << " ";
		if(r != mat.num_rows()-1) out << "| ";
	}
	out << "]";
	out << " (DenseMatrixInverse " << mat.num_rows() << "x" << mat.num_cols() << ", " << ((T::ordering == ColMajor) ? "ColMajor)" : "RowMajor)");

	return out;
}


template<typename T>
struct matrix_algebra_type_traits<DenseMatrixInverse<T> >
{
	static const int type = MATRIX_USE_GLOBAL_FUNCTIONS;
};

}  // namespace ug

#endif // not LAPACK_AVAILABLE

#endif // __H__UG__CPU_ALGEBRA__LU_DECOMP_H__
