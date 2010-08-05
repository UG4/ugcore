//	created by Andreas Vogel
//	05.08.2010

#ifndef __H__COMMON__MATH_TENSOR__
#define __H__COMMON__MATH_TENSOR__

#include <cstddef>
#include <iostream>
#include "../../types.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	MathTensor
///	a mathematical Tensor of rank TRank and N entries.

template <size_t TRank, size_t N, typename T = number> class MathTensor;

template <size_t TRank, size_t N, typename T>
class MathTensor
{
	public:
		// type of value
		typedef MathTensor<TRank-1, N, T> value_type;

		// size type
		typedef size_t size_type;

		// Dimension of tensor
		static const size_t Dimension = N;

		// Rank of tensor (i.e. number of indices)
		static const size_t Rank = TRank;

	public:
		MathTensor() {}
		MathTensor(const MathTensor& v)	{assign(v);}

		// operations with other vectors
		MathTensor& operator=  (const MathTensor& v)
		{
			if(this != &v){assign(v);}
			return *this;
		}

		inline size_t size() const {return N;}
		inline size_t rank() const {return TRank;}

		inline value_type& operator[](size_t i)				{UG_ASSERT(i < size(), "Index out of range."); return m_data[i];}
		inline const value_type& operator[](size_t i) const	{UG_ASSERT(i < size(), "Index out of range."); return m_data[i];}

	protected:
		inline void assign(const MathTensor<TRank, N, T>& v){for(size_t i = 0; i < N; ++i) m_data[i] = v[i];}

	protected:
		value_type	m_data[N];
};


template <size_t N, typename T>
class MathTensor<1, N, T>
{
	public:
		// type of value
		typedef T value_type;

		// size type
		typedef size_t size_type;

		// Dimension of tensor
		static const size_t Dimension = N;

		// Rank of tensor (i.e. number of indices)
		static const size_t Rank = 1;

	public:
		MathTensor(const bool init = true) { if(init) for(size_t i = 0; i < N; ++i) m_data[i] = 0.0;}
		MathTensor(const MathTensor& v)	{assign(v);}

		// operations with other vectors
		MathTensor& operator=  (const MathTensor& v)
		{
			if(this != &v){assign(v);}
			return *this;
		}

		inline size_t size() const {return N;}
		inline size_t rank() const {return 1;}

		inline value_type& operator[](size_t i)				{UG_ASSERT(i < size(), "Index out of range."); return m_data[i];}
		inline const value_type& operator[](size_t i) const	{UG_ASSERT(i < size(), "Index out of range."); return m_data[i];}

	protected:
		inline void assign(const MathVector<N>& v) {for(size_t i = 0; i < N; ++i) m_data[i] = v[i];}

	protected:
		value_type	m_data[N];
};

template <size_t TRank, size_t N, typename T>
std::ostream& operator<< (std::ostream& outStream, const ug::MathTensor<TRank, N, T>& v)
{
	for(size_t i = 0; i < v.size()-1; ++i)
			outStream << v[i] << ' ';
	outStream << v[v.size()-1] << ',';
	return outStream;
}

}//	end of namespace


#endif /* __H__COMMON__MATH_MathVector__ */
