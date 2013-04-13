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
/**
 * \defgroup math_tensor Tensor
 * \ingroup ugbase_math
 * \{
 */

/// a mathematical Tensor of rank TRank and N entries.
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

	///	sets all values of the tensor to a value
		MathTensor& operator=(T val) {set(val); return *this;}
		inline void set(T val) {for(size_t i = 0; i < N; ++i) m_data[i].set(val);}

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
		MathTensor(const bool init = true) { if(init) set(0.0);}
		MathTensor(const MathTensor& v)	{assign(v);}

		// operations with other vectors
		MathTensor& operator=  (const MathTensor& v)
		{
			if(this != &v){assign(v);}
			return *this;
		}

		inline size_t size() const {return N;}
		inline size_t rank() const {return 1;}

	///	sets all values of the tensor to a value
		MathTensor& operator=(T val) {set(val); return *this;}
		inline void set(T val) {for(size_t i = 0; i < N; ++i) m_data[i] = val;}

		inline value_type& operator[](size_t i)				{UG_ASSERT(i < size(), "Index out of range."); return m_data[i];}
		inline const value_type& operator[](size_t i) const	{UG_ASSERT(i < size(), "Index out of range."); return m_data[i];}

	protected:
		inline void assign(const MathTensor<1, N, T>& v) {for(size_t i = 0; i < N; ++i) m_data[i] = v[i];}

	protected:
		value_type	m_data[N];
};

template <size_t TRank, size_t N, typename T>
std::ostream& operator<< (std::ostream& outStream, const ug::MathTensor<TRank, N, T>& v)
{
	outStream << "[";
	for(size_t i = 0; i < v.size()-1; ++i)
			outStream << v[i] << ", ";
	outStream << v[v.size()-1] << "]";
	return outStream;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

// unfortunately, not all tensor are like 3x3x3x3, but also 3x2x3x2.
// merge this if possible with MathVector<MathVector<2> > ?

template <typename TEntry, size_t N>
class MathTensorX
{
public:
	// type of value
	typedef TEntry value_type;

	// size type
	typedef size_t size_type;

	// Dimension of tensor
	static const size_t Dimension = N;

	// Rank of tensor (i.e. number of indices)
	static const size_t Rank = TEntry::Rank+1;

public:
	MathTensorX() { }
	MathTensorX(const MathTensorX& v)	{assign(v);}

	// operations with other vectors
	MathTensorX& operator=  (const MathTensorX& v)
	{
		if(this != &v){assign(v);}
		return *this;
	}

	inline size_t size() const {return N;}
	inline size_t rank() const {return Rank;}

	inline value_type& operator[](size_t i)				{UG_ASSERT(i < size(), "Index out of range."); return m_data[i];}
	inline const value_type& operator[](size_t i) const	{UG_ASSERT(i < size(), "Index out of range."); return m_data[i];}

protected:
	inline void assign(const MathVector<N>& v) {for(size_t i = 0; i < N; ++i) m_data[i] = v[i];}

protected:
	TEntry m_data[N];
};

template <size_t N, typename T = number>
class MathTensor1 : public MathTensor<1, N, T> { };

template <size_t N1, size_t N2, typename T = number>
class MathTensor2 : public MathTensorX<MathTensor1<N2, T>, N1> { };

template <size_t N1, size_t N2, size_t N3, typename T = number>
class MathTensor3 : public MathTensorX< MathTensorX<MathTensor1<N3, T>, N2>, N1> { };

template <size_t N1, size_t N2, size_t N3, size_t N4, typename T = number>
class MathTensor4 : public MathTensorX< MathTensorX< MathTensorX<MathTensor1<N4, T>, N3>, N2>, N1> { };


template <typename TEntry, size_t N>
std::ostream& operator<< (std::ostream& outStream, const ug::MathTensorX<TEntry, N>& v)
{
	outStream << "[";
	for(size_t i = 0; i < N-1; ++i)
		outStream << v[i] << ", ";
	outStream << v[N-1] << "]";
	return outStream;
}

// end group math_tensor
/// \}

}//	end of namespace


#endif /* __H__COMMON__MATH_MathVector__ */
