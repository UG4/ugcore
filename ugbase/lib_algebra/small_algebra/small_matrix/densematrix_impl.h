/**
 * \file densematrix_impl.h
 *
 * \author Martin Rupp
 *
 * \date 21.07.2010
 *
 * Goethe-Center for Scientific Computing 2010.
 *
 * awesome!
 */

#ifndef __H__UG__COMMON__DENSEMATRIX_IMPL_H__
#define __H__UG__COMMON__DENSEMATRIX_IMPL_H__

#include "print.h"
#include "../blocks.h"

// constructors
namespace ug{


template<typename TStorage>
DenseMatrix<TStorage>::DenseMatrix() : TStorage()
{
}

template<typename TStorage>
DenseMatrix<TStorage>::DenseMatrix(const DenseMatrix<TStorage> &rhs) : TStorage(rhs)
{
}

template<typename TStorage>
DenseMatrix<TStorage>::DenseMatrix(double value) : TStorage()
{
	operator =(value);
}

// matrix assignment operators
////// =

template<typename TStorage>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator =(const this_type &t)
{
	resize(t.num_rows(), t.num_cols());
	for(size_t r1=0; r1<t.num_rows(); r1++)
		for(size_t c1=0; c1<t.num_cols(); c1++)
			entry(r1, c1) = t(r1, c1);
	return *this;
}

template<typename TStorage>
template<typename T2>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator =(const T2 &t)
{
	resize(t.num_rows(), t.num_cols());
	for(size_t r1=0; r1<t.num_rows(); r1++)
		for(size_t c1=0; c1<t.num_cols(); c1++)
			entry(r1, c1) = t(r1, c1);
	return *this;
}

template<typename TStorage>
DenseMatrix<TStorage> &
// this is stupid since i like to make it that way that i have one ::value_type and one double,
// but then there are two (double) functions...
//DenseMatrix<TStorage>::operator = (const typename DenseMatrix<TStorage>::value_type &rhs)
DenseMatrix<TStorage>::operator = (double rhs)
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
		{
			if(r==c) entry(r, c) = rhs;
			else     entry(r, c) = 0.0;
		}
	return *this;
}

////// +=
template<typename TStorage>
template<typename T2>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator += (const T2 &t)
{
	UG_ASSERT(t.num_rows() == num_rows() && t.num_cols() == num_cols(), "");
	for(size_t r1=0; r1<t.num_rows(); r1++)
		for(size_t c1=0; c1<t.num_cols(); c1++)
			entry(r1, c1) += t(r1, c1);
	return *this;
}

template<typename TStorage>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator += (double alpha)
{
	size_t minimum=num_rows() > num_cols() ? num_cols() : num_rows();
	for(size_t i=0; i<minimum; i++)
			entry(i, i) += alpha;
	return *this;
}

////// -=

template<typename TStorage>
template<typename T2>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator -= (const T2 &t)
{
	UG_ASSERT(t.num_rows() == num_rows() && t.num_cols() == num_cols(), "");
	for(size_t r1=0; r1<t.num_rows(); r1++)
		for(size_t c1=0; c1<t.num_cols(); c1++)
			entry(r1, c1) -= t(r1, c1);
	return *this;
}

template<typename TStorage>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator -= (double alpha)
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			entry(r, c) -= alpha;
	return *this;
}

////// *=

template<typename TStorage>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator *= (double alpha)
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			entry(r, c) *= alpha;
	return *this;
}

template<typename TStorage>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator *= (const this_type &mat)
{
	operator=(operator*(mat));
	return *this;
}

////// /=
template<typename TStorage>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator /= (double alpha)
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			entry(r, c) /= alpha;
	return *this;
}

template<typename TStorage>
DenseMatrix<TStorage> &
DenseMatrix<TStorage>::operator /= (this_type &other)
{
	this_type tmp = other;
	Invert(tmp);
	(*this) = (*this) * tmp;
	return *this;
}


////// +
template<typename TStorage>
DenseMatrix<TStorage> 
DenseMatrix<TStorage>::operator + (const this_type &other ) const
{
	UG_ASSERT(num_rows() == other.num_rows() && num_cols() == other.num_cols(), "");
	this_type erg;
	erg.resize(num_rows(), num_cols());
	for(size_t r=0; r<num_rows(); r++)
		for(size_t c=0; c<num_cols(); c++)
			erg(r, c) = entry(r, c) + other(r,c);
	return erg;
}
////// -
template<typename TStorage>
DenseMatrix<TStorage> 
DenseMatrix<TStorage>::operator - (const this_type &other ) const
{
	UG_ASSERT(num_rows() == other.num_rows() && num_cols() == other.num_cols(), "");
	this_type erg;
	erg.resize(num_rows(), num_cols());

	for(size_t r=0; r<num_rows(); r++)
		for(size_t c=0; c<num_cols(); c++)
			erg(r, c) = entry(r, c) - other(r,c);
	return erg;
}

////// unary -
template<typename TStorage>
DenseMatrix<TStorage> 
DenseMatrix<TStorage>::operator - () const
{
	this_type erg;
	erg.resize(num_rows(), num_cols());
	for(size_t r=0; r < num_rows(); r++)
		for(size_t c=0; c < num_cols(); c++)
		{
			erg(r,c) = entry(r, c);
			erg(r,c) *= -1.0;
		}
	return erg;
}

// multiply
////// *
template<typename TStorage>
DenseMatrix<TStorage> 
DenseMatrix<TStorage>::operator * (const this_type &other ) const
{
	// that aint 100% correct
	UG_ASSERT(num_cols() == other.num_rows(), "");

	this_type erg;
	erg.resize(num_rows(), other.num_cols());

	for(size_t r=0; r < num_rows(); r++)
		for(size_t c=0; c < other.num_cols(); c++)
		{
			erg(r,c) = 0.0;
			for(size_t i=0; i < num_cols(); i++)
				AddMult(erg(r,c), at(r, i), other.at(i, c));
		}
	return erg;
}


template<typename TStorage>
DenseMatrix<TStorage>
DenseMatrix<TStorage>::T() const
{
	this_type erg;
	erg.resize(num_rows(), num_cols());
	for(size_t r=0; r < num_rows(); r++)
		for(size_t c=0; c < num_cols(); c++)
			erg(r,c) = entry(c, r);
	return erg;
}

template<typename TStorage>
template<typename TStorage2>
DenseVector<TStorage2>
DenseMatrix<TStorage>::operator * (const DenseVector<TStorage2> &vec) const
{
	UG_ASSERT(num_cols() == vec.size(), "");
	DenseVector<TStorage2> erg;
	erg.resize(num_rows());

	for(size_t r=0; r < num_rows(); r++)
	{
		erg[r] = 0.0;
		for(size_t c=0; c < num_cols(); c++)
			erg[r] += at(r,c) * vec[c];
	}
	return erg;
}

template<typename TStorage>
DenseMatrix<TStorage> 
DenseMatrix<TStorage>::operator * (double alpha ) const
{
	this_type erg;
	erg.resize(num_rows(), num_cols());

	for(size_t r=0; r < num_rows(); r++)
		for(size_t c=0; c < num_cols(); c++)
			erg(r,c) = at(r,c)*alpha;
	return erg;
}



///// /
template<typename TStorage>
DenseMatrix<TStorage> 
DenseMatrix<TStorage>::operator / (this_type &other)
{
	this_type tmp = other;
	Invert(tmp);

	return (*this) * tmp;
}

// compare operators

template<typename TStorage>
bool
DenseMatrix<TStorage>::operator == (double t) const
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
		{
			if(r==c)
			{
				if(entry(r,c) != t) return false;
			}
			else
				if(entry(r,c) != 0.0) return false;
		}
	return true;
}

template<typename TStorage>
template<typename TStorage2>
bool
DenseMatrix<TStorage>::operator == (const DenseMatrix<TStorage2> &t) const
{
	for(size_t r=0; r<num_rows(); ++r)
		for(size_t c=0; c<num_cols(); ++c)
			if(entry(r,c) != t(r,c)) return false;
	return true;
}


template<typename TStorage>
template<typename T2>
bool DenseMatrix<TStorage>::operator != (const T2 &t) const
{
	return !(operator == (t));
}

template<typename TStorage>
void DenseMatrix<TStorage>::maple_print(const char *name)
{
	UG_LOG(MatlabString(*this, name));
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
//	out << "(DenseMatrix " << mat.num_rows() << "x" << mat.num_cols() << ", " << ((DenseMatrix<TStorage>::ordering == ColMajor) ? "ColMajor)" : "RowMajor)");

	return out;
}


template<size_t Tr, size_t Tc>
inline
bool
BlockSerialize(const DenseMatrix<FixedArray2<number, Tr, Tc> > &mat, std::ostream &buff)
{
	buff.write((char*)&mat, sizeof(mat));
	return true;
}

template<size_t Tr, size_t Tc>
inline
bool
BlockDeserialize(std::istream &buff, const DenseMatrix<FixedArray2<number, Tr, Tc> > &mat)
{
	buff.read((char*)&mat, sizeof(mat));
	return true;
}


template<typename T>
inline
void
Serialize(std::ostream &buff, const DenseMatrix<VariableArray2<T> > &mat)
{
	size_t rows = mat.num_rows();
	size_t cols = mat.num_cols();
	buff.write((char*)&rows, sizeof(rows));
	buff.write((char*)&cols, sizeof(cols));
	for(size_t r=0; r<rows; r++)
		for(size_t c=0; c<cols; c++)
			BlockSerialize(mat(r, c), buff);
}


template<typename T>
inline
void
Deserialize(std::istream &buff, const DenseMatrix<VariableArray2<T> > &mat)
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
