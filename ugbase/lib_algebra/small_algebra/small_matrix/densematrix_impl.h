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
	if(this == &rhs) return *this;

	if(num_rows() != rhs.num_rows() || num_cols() != rhs.num_cols())
		resize(rhs.num_rows(), rhs.num_cols());

	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			entry(r, c) = rhs(r, c);
	return *this;
}

template<typename TStorage>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator += (const DenseMatrix<TStorage> &rhs)
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			entry(r, c) += rhs(r, c);
	return *this;
}

template<typename TStorage>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator -= (const DenseMatrix<TStorage> &rhs)
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			entry(r, c) -= rhs(r, c);
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
			entry(r, c) = (r==c ? rhs : 0.0);
	return *this;
}

template<typename TStorage>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator += (const DenseMatrix<TStorage>::value_type &alpha)
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			entry(r, c) += alpha;
	return *this;
}

template<typename TStorage>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator -= (const DenseMatrix<TStorage>::value_type &alpha)
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			entry(r, c) -= alpha;
	return *this;
}

template<typename TStorage>
template<typename T>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator *= (const T &alpha)
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			entry(r, c) *= alpha;
	return *this;
}

template<typename TStorage>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator /= (const DenseMatrix<TStorage>::value_type &alpha)
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			entry(r, c) /= alpha;
	return *this;
}

// compare operators

template<typename TStorage>
template<typename T>
bool DenseMatrix<TStorage>::operator == (const T &t) const
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			if(entry(r,c) != t) return false;
	return true;
}

template<typename TStorage>
template<typename T>
bool DenseMatrix<TStorage>::operator == (const DenseMatrix<T> &t) const
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			if(entry(r,c) != t(r,c)) return false;
	return true;
}


template<typename TStorage>
template<typename T>
bool DenseMatrix<TStorage>::operator != (const T &t) const
{
	return !(operator == (t));
}

template<typename TStorage>
void DenseMatrix<TStorage>::maple_print(const char *name)
{
	UG_LOG(name << " = matrix([");
	for(size_t r=0; r<num_rows(); ++r)
	{
		if(r > 0) UG_LOG(", ");
		UG_LOG("[");
		for(size_t c=0; c<num_cols(); ++c)
		{
			if(c > 0) UG_LOG(", ");
			UG_LOG(entry(r, c));
		}
		UG_LOG("]");
	}
	UG_LOG("]);\n");
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


template<size_t Tr, size_t Tc> inline bool BlockSerialize(const DenseMatrix<FixedArray2<number, Tr, Tc> > &mat, std::ostream &buff)
{
	buff.write((char*)&mat, sizeof(mat));
	return true;
}

template<size_t Tr, size_t Tc> inline bool BlockDeserialize(std::istream &buff, const DenseMatrix<FixedArray2<number, Tr, Tc> > &mat)
{
	buff.read((char*)&mat, sizeof(mat));
	return true;
}


template<typename T> inline void Serialize(std::ostream &buff, const DenseMatrix<VariableArray2<T> > &mat)
{
	size_t rows = mat.num_rows();
	size_t cols = mat.num_cols();
	buff.write((char*)&rows, sizeof(rows));
	buff.write((char*)&cols, sizeof(cols));
	for(size_t r=0; r<rows; r++)
		for(size_t c=0; c<cols; c++)
			BlockSerialize(mat(r, c), buff);
}

template<typename T> inline void Deserialize(std::istream &buff, const DenseMatrix<VariableArray2<T> > &mat)
{
	size_t rows, cols;
	buff.read((char*)&rows, sizeof(rows));
	buff.read((char*)&cols, sizeof(cols));
	mat.resize(rows, cols);
	for(size_t r=0; r<rows; r++)
		for(size_t c=0; c<cols; c++)
			BlockDeserialize(buff, mat(r, c));
}


} // namespace ug

#endif // __H__UG__COMMON__DENSEMATRIX_IMPL_H__
