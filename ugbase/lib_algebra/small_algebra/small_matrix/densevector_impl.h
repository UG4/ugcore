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
		entry(i) = 0.0;
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
	if(this == &rhs) return *this;
	if(size() != rhs.size()) resize(rhs.size());

	for(size_type i=0; i<size(); i++)
		entry(i) = rhs[i];
	return *this;
}


template<typename TStorage>
DenseVector<TStorage> &
DenseVector<TStorage>::operator += (const DenseVector<TStorage> &rhs)
{
	for(size_type i=0; i<size(); i++)
		entry(i) += rhs[i];
	return *this;
}


template<typename TStorage>
DenseVector<TStorage> &
DenseVector<TStorage>::operator -= (const DenseVector<TStorage> &rhs)
{
	for(size_type i=0; i<size(); i++)
		entry(i) -= rhs[i];
	return *this;
}


// operations with scalars
template<typename TStorage>
template<typename T>
DenseVector<TStorage> &
DenseVector<TStorage>::operator=(const T &alpha)
{
	for(size_t i=0; i<size(); i++)
		entry(i) = alpha;
	return *this;
}

template<typename TStorage>
DenseVector<TStorage> &
DenseVector<TStorage>::operator+=(const typename DenseVector<TStorage>::value_type &alpha)
{
	for(size_t i=0; i<size(); i++)
		entry(i) += alpha;
	return *this;
}

template<typename TStorage>
DenseVector<TStorage> &
DenseVector<TStorage>::operator-=(const typename DenseVector<TStorage>::value_type &alpha)
{
	for(size_t i=0; i<size(); i++)
		entry(i) -= alpha;
	return *this;
}

template<typename TStorage>
template<typename T>
DenseVector<TStorage> &
DenseVector<TStorage>::operator*=(const T &alpha)
{
	for(size_t i=0; i<size(); i++)
		entry(i) *= alpha;
	return *this;
}


template<typename TStorage>
DenseVector<TStorage> &
DenseVector<TStorage>::operator/=(const typename DenseVector<TStorage>::value_type &alpha)
{
	for(size_t i=0; i<size(); i++)
		entry(i) /= alpha;

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



template<typename TStorage>
template<typename Type>
DenseVector<TStorage> &
DenseVector<TStorage>::assign(const Type &t)
{
	VectorAssign(*this, t);
	return *this;
}



template<size_t n> inline bool BlockSerialize(const DenseVector<FixedArray1<number, n> > &vec, std::ostream &buff)
{
	buff.write((char*)&vec, sizeof(vec));
	return true;
}

template<size_t n> inline bool BlockDeserialize(std::istream &buff, DenseVector<FixedArray1<number, n> > &vec)
{
	buff.read((char*)&vec, sizeof(vec));
	return true;
}


template<typename T> inline bool BlockSerialize(const DenseVector<VariableArray1<T> > &vec, std::ostream &buff)
{
	size_t s = vec.size();
	buff.write((char*)&s, sizeof(s));
	for(size_t i=0; i<s; i++)
		BlockSerialize(vec[i], buff);
	return true;
}

template<typename T> inline bool BlockDeserialize(std::istream &buff, DenseVector<VariableArray1<double> > &vec)
{
	size_t s;
	buff.read((char*)&s, sizeof(s));
	vec.resize(s);
	for(size_t i=0; i<s; i++)
		BlockDeserialize(buff, vec[i]);
	return true;
}

//MAKE_TEMPLATE_OPERATORS_VECTOR2(typename TStorage, DenseVector<TStorage>);

}

#endif // __H__UG__COMMON__DENSEVECTOR_IMPL_H__
