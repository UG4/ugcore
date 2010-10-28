/**
 * \file densevector_impl.h
 *
 * \author Martin Rupp
 *
 * \date 21.07.2010
 *
 * Goethe-Center for Scientific Computing 2010.
 */


#ifndef __H__UG__COMMON__DENSEVECTOR_IMPL_H__
#define __H__UG__COMMON__DENSEVECTOR_IMPL_H__

namespace ug{

template<typename TStorage>
DenseVector<TStorage>::DenseVector()
{
	for(size_t i=0; i<size(); i++)
		at(i) = 0.0;
}

template<typename TStorage>
DenseVector<TStorage>::DenseVector(const DenseVector<TStorage> &rhs) : TStorage(rhs)
{
}

// operations with vectors
template<typename TStorage>
DenseVector<TStorage> &
DenseVector<TStorage>::operator = (const DenseVector<TStorage> &rhs)
{
	for(size_type i=0; i<size(); i++)
		at(i) = rhs[i];
	return *this;
}


template<typename TStorage>
DenseVector<TStorage> &
DenseVector<TStorage>::operator += (const DenseVector<TStorage> &rhs)
{
	for(size_type i=0; i<size(); i++)
		at(i) += rhs[i];
	return *this;
}


template<typename TStorage>
DenseVector<TStorage> &
DenseVector<TStorage>::operator -= (const DenseVector<TStorage> &rhs)
{
	for(size_type i=0; i<size(); i++)
		at(i) -= rhs[i];
	return *this;
}


// operations with scalars
template<typename TStorage>
template<typename T>
DenseVector<TStorage> &
DenseVector<TStorage>::operator=(const T &alpha)
{
	for(size_t i=0; i<size(); i++)
		at(i) = alpha;
	return *this;
}

template<typename TStorage>
DenseVector<TStorage> &
DenseVector<TStorage>::operator+=(const typename DenseVector<TStorage>::value_type &alpha)
{
	for(size_t i=0; i<size(); i++)
		at(i) += alpha;
	return *this;
}

template<typename TStorage>
DenseVector<TStorage> &
DenseVector<TStorage>::operator-=(const typename DenseVector<TStorage>::value_type &alpha)
{
	for(size_t i=0; i<size(); i++)
		at(i) -= alpha;
	return *this;
}

template<typename TStorage>
template<typename T>
DenseVector<TStorage> &
DenseVector<TStorage>::operator*=(const T &alpha)
{
	for(size_t i=0; i<size(); i++)
		at(i) *= alpha;
	return *this;
}


template<typename TStorage>
DenseVector<TStorage> &
DenseVector<TStorage>::operator/=(const typename DenseVector<TStorage>::value_type &alpha)
{
	for(size_t i=0; i<size(); i++)
		at(i) /= alpha;

	return *this;
}

// views
// methods


template<typename TStorage>
std::ostream &operator << (std::ostream &out, const DenseVector<TStorage> &vec)
{
	out << "[";
	for(size_t i=0; i<vec.size(); i++)
		out << " " << vec[i];
	out << " ] (" << vec.size() << ")";
	return out;
}

}

#endif // __H__UG__COMMON__DENSEVECTOR_IMPL_H__
