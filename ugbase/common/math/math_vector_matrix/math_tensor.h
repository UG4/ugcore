/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

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
	inline void assign(const MathTensorX<TEntry, N>& v) {for(size_t i = 0; i < N; ++i) m_data[i] = v[i];}

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

template <typename TEntry>
std::ostream& operator<< (std::ostream& outStream, const ug::MathTensorX<TEntry, 1>& v)
{
	outStream << "[";
	outStream << v[0] << "]";
	return outStream;
}

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
