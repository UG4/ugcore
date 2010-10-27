#ifndef __H__UG__CPU_ALGEBRA__LAPACK_INVERT_H__
#define __H__UG__CPU_ALGEBRA__LAPACK_INVERT_H__

namespace ug{

//////////////////////////////////////////////////////
// 1x1

template<typename T>
inline bool GetInverse1(DenseMatrix<T> &inv, const DenseMatrix<T> &mat)
{
	UG_ASSERT(mat(0,0)!=0.0, "Determinant zero, cannot invert matrix.");
	if(mat(0,0) == 0.0) return false;
	inv(0,0) = 1/mat(0,0);
	return true;
}

template<typename T>
bool Invert1(DenseMatrix<T> &mat)
{
	UG_ASSERT(mat(0,0)!=0.0, "Determinant zero, cannot invert matrix.");
	if(mat(0,0) == 0.0) return false;
	mat(0,0) = 1/mat(0,0);
	return true;
};

inline bool GetInverse(DenseMatrix<FixedArray2<double, 1, 1> > &inv, const DenseMatrix<FixedArray2<double, 1, 1> > &mat)
{
	return GetInverse1(inv, mat);
}

inline bool Invert(DenseMatrix< FixedArray2<double, 1, 1> > &mat)
{
	return Invert1(mat);
}


//////////////////////
// 2x2

template<typename T>
inline double GetDet2(const DenseMatrix<T> &mat)
{
	UG_ASSERT(mat.num_rows() == 2 && mat.num_cols() == 2, "only for 2x2-matrices");
	return mat(0,0)*mat(1,1) - mat(1,0)*mat(0,1);
}

template<typename T>
inline bool GetInverse2(DenseMatrix<T> &inv, const DenseMatrix<T> &mat)
{
	UG_ASSERT(&inv != &mat, "inv and mat have to be different. Otherwise use Invert/Invert2");
	double invdet = GetDet2(mat);
	UG_ASSERT(invdet != 0, "Determinant zero, cannot invert matrix.");
	if(invdet == 0.0) return false;
	invdet = 1.0/invdet;
	inv(0,0) = mat(1,1) * invdet;
	inv(1,1) = mat(0,0) * -invdet;
	inv(0,1) = mat(1,0) * -invdet;
	inv(1,0) = mat(0,1) * invdet;
	return true;
}

template<typename T>
bool Invert2(DenseMatrix<T> &mat)
{

	double invdet = GetDet2(mat);
	UG_ASSERT(invdet != 0, "Determinant zero, cannot invert matrix.");
	if(invdet == 0.0) return false;
	invdet = 1.0/invdet;

	swap(mat(0,0), mat(1,1));
	swap(mat(1,0), mat(0,1));

	mat(0,0) *= invdet;
	mat(0,1) *= -invdet;
	mat(1,0) *= -invdet;
	mat(1,1) *= invdet;
	return true;
};

template<>
inline bool GetInverse(DenseMatrix<FixedArray2<double, 2, 2> > &inv, const DenseMatrix<FixedArray2<double, 2, 2> > &mat)
{
	return GetInverse2(inv, mat);
}

template<>
inline bool Invert(DenseMatrix< FixedArray2<double, 2, 2> > &mat)
{
	return Invert2(mat);
}


//////////////////////
// 3x3


template<typename T>
inline double GetDet3(const DenseMatrix<T> &mat)
{
	UG_ASSERT(mat.num_rows() == 3 && mat.num_cols() == 3, "only for 3x3-matrices");
	return 	mat(0,0)*mat(1,1)*mat(2,2) + mat(0,1)*mat(1,2)*mat(2,0) + mat(0,2)*mat(1,0)*mat(2,1)
			- mat(0,0)*mat(1,2)*mat(2,1) - mat(0,1)*mat(1,0)*mat(2,2) - mat(0,2)*mat(1,1)*mat(2,0);
}

template<typename T>
inline bool GetInverse3(DenseMatrix<T> &inv, const DenseMatrix<T> &mat)
{
	UG_ASSERT(&inv != &mat, "inv and mat have to be different. Otherwise use Invert/Invert3");
	double invdet = GetDet3(mat);
	UG_ASSERT(invdet != 0, "Determinant zero, cannot invert matrix.");
	if(invdet == 0.0) return false;
	invdet = 1.0/invdet;

	inv(0,0) = ( mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)) * invdet;
	inv(0,1) = (-mat(0,1)*mat(2,2) + mat(0,2)*mat(2,1)) * invdet;
	inv(0,2) = ( mat(0,1)*mat(1,2) - mat(0,2)*mat(1,1)) * invdet;
	inv(1,0) = (-mat(1,0)*mat(2,2) + mat(1,2)*mat(2,0)) * invdet;
	inv(1,1) = ( mat(0,0)*mat(2,2) - mat(0,2)*mat(2,0)) * invdet;
	inv(1,2) = (-mat(0,0)*mat(1,2) + mat(0,2)*mat(1,0)) * invdet;
	inv(2,0) = ( mat(1,0)*mat(2,1) - mat(1,1)*mat(2,0)) * invdet;
	inv(2,1) = (-mat(0,0)*mat(2,1) + mat(0,1)*mat(2,0)) * invdet;
	inv(2,2) = ( mat(0,0)*mat(1,1) - mat(0,1)*mat(1,0)) * invdet;
	return true;
}

inline bool Invert3(DenseMatrix<FixedArray2<double, 3, 3> > & mat)
{
	DenseMatrix<FixedArray2<double, 3, 3> > inv;
	if(GetInverse3(inv, mat) == false) return false;
	mat = inv;
	return true;
}

inline bool Invert3(DenseMatrix<VariableArray2<double> > & mat)
{
	DenseMatrix<VariableArray2<double> > inv;
	inv.resize(3,3);
	if(GetInverse3(inv, mat) == false) return false;
	mat = inv;
	return true;
}

template<>
inline bool GetInverse(DenseMatrix<FixedArray2<double, 3, 3> > &inv, const DenseMatrix<FixedArray2<double, 3, 3> > &mat)
{
	return GetInverse3(inv, mat);
}

template<>
inline bool Invert(DenseMatrix< FixedArray2<double, 3, 3> > &mat)
{
	return Invert3(mat);
}

//////////////////////

template<typename T>
void InvertN(DenseMatrix<T> &mat, static_type /*dummy*/)
{

}

template<typename T>
bool InvertNdyn(DenseMatrix<T> &mat)
{
	std::vector<lapack_int> interchange;
	interchange.resize(mat.num_rows());

	int info = getrf(mat.num_rows(), mat.num_cols(), &mat(0,0), mat.num_rows(), &interchange[0]);
	UG_ASSERT(info == 0, "info is " << info << ( info > 0 ? ": SparseMatrix singular in U(i,i)" : ": i-th argument had had illegal value"));
	if(info == 0) return false;

	// calc work size
	double worksize; int iWorksize = -1;
	info = getri(mat.num_rows(), &mat(0,0), mat.num_rows(), &interchange[0], &worksize, iWorksize);
	UG_ASSERT(info == 0, "");
	iWorksize = worksize;

	std::vector<double> work;
	work.resize(iWorksize);

	info = getri(mat.num_rows(), &mat(0,0), mat.num_rows(), &interchange[0], &work[0], iWorksize);
	UG_ASSERT(info == 0, "");
	if(info == 0) return false;

	return true;
}

template<typename T, size_t TUnknowns>
bool Invert(DenseMatrix<FixedArray2<T, TUnknowns, TUnknowns> > &mat)
{
	lapack_int interchange[T::static_row_size];

	int info = getrf(mat.num_rows(), mat.num_cols(), &mat(0,0), mat.num_rows(), interchange);
	UG_ASSERT(info == 0, "info is " << info << ( info > 0 ? ": Matrix singular in mat(i,i)" : ": i-th argument had had illegal value"));
	if(info == 0) return false;

	// calc work size
	// todo: make this static
	double worksize; int iWorksize = -1;
	info = getri(mat.num_rows(), &mat(0,0), mat.num_rows(), interchange, &worksize, iWorksize);
	UG_ASSERT(info == 0, "");
	iWorksize = worksize;

	std::vector<double> work;
	work.resize(iWorksize);

	info = getri(mat.num_rows(), mat(0,0), mat.num_rows(), interchange, &work[0], iWorksize);
	UG_ASSERT(info == 0, "");
	if(info == 0) return false;

	return true;
}

template<typename T>
inline bool Invert(DenseMatrix<T> &mat)
{
	switch(mat.num_rows())
	{
		case 1: return Invert1(mat);
		case 2: return Invert2(mat);
		case 3: return Invert3(mat);
		default: return InvertNdyn(mat);
	}
}

}
#endif
