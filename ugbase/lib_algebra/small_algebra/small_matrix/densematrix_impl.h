/**
 * \file densematrix_impl.h
 *
 * \author Martin Rupp
 *
 * \date 21.07.2010
 *
 * Goethe-Center for Scientific Computing 2010.
 */

#ifndef __H__UG__COMMON__DENSEMATRIX_IMPL_H__
#define __H__UG__COMMON__DENSEMATRIX_IMPL_H__

// constructors
namespace ug{


template<typename TStorage>
DenseMatrix<TStorage>::DenseMatrix()
{
}

template<typename TStorage>
DenseMatrix<TStorage>::DenseMatrix(const DenseMatrix<TStorage> &rhs) : TStorage(rhs)
{
}

// matrix assignment operators

template<typename TStorage>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator = (const DenseMatrix<TStorage> &rhs)
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			at(r, c) = rhs(r, c);
	return *this;
}

template<typename TStorage>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator += (const DenseMatrix<TStorage> &rhs)
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			at(r, c) += rhs(r, c);
	return *this;
}

template<typename TStorage>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator -= (const DenseMatrix<TStorage> &rhs)
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			at(r, c) -= rhs(r, c);
	return *this;
}

// alpha operators

template<typename TStorage>
template<typename T>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator = (const T &rhs)
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			at(r, c) = (r==c ? rhs : 0.0);
	return *this;
}

template<typename TStorage>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator += (const DenseMatrix<TStorage>::value_type &alpha)
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			at(r, c) += alpha;
	return *this;
}

template<typename TStorage>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator -= (const DenseMatrix<TStorage>::value_type &alpha)
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			at(r, c) -= alpha;
	return *this;
}

template<typename TStorage>
template<typename T>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator *= (const T &alpha)
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			at(r, c) *= alpha;
	return *this;
}

template<typename TStorage>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator /= (const DenseMatrix<TStorage>::value_type &alpha)
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			at(r, c) /= alpha;
	return *this;
}

// compare operators

template<typename TStorage>
template<typename T>
bool DenseMatrix<TStorage>::operator == (const T &t) const
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			if(at(r,c) != t) return false;
	return true;
}

template<typename TStorage>
template<typename T>
bool DenseMatrix<TStorage>::operator == (const DenseMatrix<T> &t) const
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			if(at(r,c) != t(r,c)) return false;
	return true;
}


template<typename TStorage>
template<typename T>
bool DenseMatrix<TStorage>::operator != (const T &t) const
{
	return !(operator == (t));
}


template<typename TStorage>
std::ostream &operator << (std::ostream &out, const DenseMatrix<TStorage> &mat)
{
	out << "[ ";
	typedef size_t size_type;
	for(size_type r=0; r<mat.num_rows(); ++r)
	{
		for(size_type c=0; c<mat.num_cols(); ++c)
			out << mat(r, c) << " ";
		if(r != mat.num_rows()-1) out << "| ";
	}
	out << "]";
	out << "(DenseMatrix " << mat.num_rows() << "x" << mat.num_cols() << ", " << ((DenseMatrix<TStorage>::ordering == ColMajor) ? "ColMajor)" : "RowMajor)");

	return out;
}



} // namespace ug

#endif // __H__UG__COMMON__DENSEMATRIX_IMPL_H__
