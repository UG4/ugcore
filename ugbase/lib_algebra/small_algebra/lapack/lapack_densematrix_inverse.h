#ifndef __H__UG__SMALL_ALGEBRA__DENSEMATRIX_INVERT_H__
#define __H__UG__SMALL_ALGEBRA__DENSEMATRIX_INVERT_H__

#include "../small_matrix/densematrix.h"
#include "../small_matrix/densevector.h"
#include "../../common/operations.h"

namespace ug
{
///////////////////////////////////////////////////////////////////////////////////////
/**
 * smallInverse<size_t n>
 * A class to hold a inverse of a smallMatrix<n>
 * implemented with LAPACKs LU-Decomposition dgetrf
 * (uses double[n*n] for LU and interchange[n] for pivoting
 */
template<typename TStorage>
class DenseMatrixInverse
{
private: // storage
	DenseMatrix<TStorage> densemat;
	std::vector<lapack_int> interchange;

public:
	inline size_t num_cols() const
	{
		return densemat.num_cols();
	}

	inline size_t num_rows() const
	{
		return densemat.num_rows();
	}

	inline void resize(size_t k)
	{
		densemat.resize(k,k);
		densemat = 0.0;
	}
	inline void resize(size_t k, size_t k2)
	{
		UG_COND_THROW(k!=k2, "only square matrices supported");
		resize(k);
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
		if(densemat.num_rows() == 0) return false;

		interchange.resize(densemat.num_rows());
		int info = getrf(densemat.num_rows(), densemat.num_cols(), &densemat(0,0),
				densemat.num_rows(), &interchange[0]);
		if(info != 0)
		{
			if(info > 0)
			{UG_THROW("ERROR in 'DenseMatrixInverse::invert': Matrix singular in U(i,i), with i="<<info<<"\n");}
			if(info < 0)
			{UG_THROW("ERROR in 'DenseMatrixInverse::invert':  i-th argument had had illegal value, with i="<<info<<"\n");}
		}
		return info == 0;
	}

	template<typename vector_t>
	void apply(DenseVector<vector_t> &vec) const
	{
		if(vec.size() == 0) return;
		int info = getrs(ModeNoTrans, num_rows(), 1, &densemat(0,0), num_rows(), &interchange[0], &vec[0], num_rows());
		(void) info;
		UG_ASSERT(info == 0, "DenseMatrixInverse::mat_mult: getrs failed.");
	}

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

}
#endif
