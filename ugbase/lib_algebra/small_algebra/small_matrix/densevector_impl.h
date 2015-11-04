

#ifndef __H__UG__COMMON__DENSEVECTOR_IMPL_H__
#define __H__UG__COMMON__DENSEVECTOR_IMPL_H__

#include "common/serialization.h"

namespace ug{

template<typename TStorage>
DenseVector<TStorage>::DenseVector(double alpha)
{
	operator=(alpha);
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

template<typename TStorage>
void DenseVector<TStorage>::maple_print(const char *name)
{
	UG_LOG(name << " = vector([");
	for(size_t i=0; i<size(); ++i)
	{
		if(i > 0) UG_LOG(", ");
		UG_LOG(entry(i));
	}
	UG_LOG("]);\n");
}

// views
// methods


template<typename TStorage>
std::ostream &operator << (std::ostream &out, const DenseVector<TStorage> &vec)
{
	out << "[";
	for(size_t i=0; i<vec.size(); i++)
		out << " " << vec[i];
	out << " ] ";
//	out << "(" << vec.size() << ")";
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



template<size_t n> inline void Serialize(std::ostream &buff, const DenseVector<FixedArray1<number, n> > &vec)
{
	buff.write((char*)&vec, sizeof(vec));
}

template<size_t n> inline void Deserialize(std::istream &buff, DenseVector<FixedArray1<number, n> > &vec)
{
	buff.read((char*)&vec, sizeof(vec));
}


template<typename T> inline void Serialize(std::ostream &buff, const DenseVector<VariableArray1<T> > &vec)
{
	size_t s = vec.size();
	buff.write((char*)&s, sizeof(s));
	for(size_t i=0; i<s; i++)
		Serialize(buff, vec[i]);
}

template<typename T> inline void Deserialize(std::istream &buff, DenseVector<VariableArray1<double> > &vec)
{
	size_t s;
	buff.read((char*)&s, sizeof(s));
	vec.resize(s);
	for(size_t i=0; i<s; i++)
		Deserialize(buff, vec[i]);
}


template<typename T >
inline bool IsFiniteAndNotTooBig(const DenseVector<T> &v)
{
	for(size_t r=0; r<v.size(); r++)
		if(IsFiniteAndNotTooBig(v[r]) == false) return false;

	return true;
}


//MAKE_TEMPLATE_OPERATORS_VECTOR2(typename TStorage, DenseVector<TStorage>);

}

#endif // __H__UG__COMMON__DENSEVECTOR_IMPL_H__
