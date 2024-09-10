/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	This header defines common vector-types.
//	It is possible to completely avoid these vectors and to use your own.
//	Have a look at lgmath.h to see which typedefs have to be replaced.
//	You have to make sure that your vector-types specialize the
//	template methods defined in lgmath_vector_descriptor.

#ifndef __H__COMMON__MATH_VECTOR__
#define __H__COMMON__MATH_VECTOR__

#include <cstddef>
#include <iostream>
#include <algorithm>
#include "../../ug_config.h"
#include "../../types.h"
#include "common/math/misc/math_constants.h"


namespace ug
{

/**
 * \defgroup vectors Vectors
 * Abbreviations of small vectors
 * \ingroup ugbase_math
 * \{
 */

////////////////////////////////////////////////////////////////////////
//	MathMathVector
///	a mathematical Vector with N entries.
template <std::size_t N, typename T = number> class MathVector;

///	helper method which creates a vector from another vector of different dimensionality
/** Typically, the method isn't invoked directly but serves as an implementation
 * helper for MathVector<toN, T>::from(v).
 * \{ */
template <std::size_t fromN, std::size_t toN, typename T>
MathVector<toN, T> MathVectorFrom (const MathVector<fromN, T>& v)
{
	MathVector<toN, T> r;
	static const size_t minN = std::min(toN, fromN);
	static const size_t maxN = std::max(toN, fromN);
	for(size_t i = 0; i < minN; ++i)
		r[i] = v[i];
	for(size_t i = minN; i < maxN; ++i)
		r[i] = 0;
	return r;
}

template <std::size_t N, typename T>
MathVector<N, T> MathVectorFrom (const MathVector<N, T>& v)
{
	return v;
}
/** \} */


/**
 * A mathematical Vector with N entries and static storage
 */
template <std::size_t N, typename T>
class MathVector
{
	public:
		typedef T value_type;
		typedef std::size_t size_type;
		static const std::size_t Size = N;

	public:
		MathVector() {for(std::size_t i = 0; i < N; ++i) m_data[i] = 0.0;}
		MathVector(const value_type& val) {for(std::size_t i = 0; i < N; ++i) m_data[i] =  val;}
		MathVector(const MathVector& v)	{assign(v);}

		template <std::size_t fromN>
		static inline MathVector from(const MathVector<fromN, T>& v)
		{
			return MathVectorFrom<fromN, N, T>(v);
		}

		// operations with other vectors
		MathVector& operator=  (const MathVector& v)
		{
			if(this != &v)
			{
				assign(v);
			}
			return *this;
		}
		MathVector& operator+= (const MathVector& v) {for(std::size_t i = 0; i < N; ++i) m_data[i] += v.coord(i);return *this;}
		MathVector& operator-= (const MathVector& v) {for(std::size_t i = 0; i < N; ++i) m_data[i] -= v.coord(i);return *this;}

		// operations with scalar
		MathVector& operator=  (const value_type& val) {for(std::size_t i = 0; i < N; ++i) m_data[i] =  val;return *this;}
		MathVector& operator+= (const value_type& val) {for(std::size_t i = 0; i < N; ++i) m_data[i] += val;return *this;}
		MathVector& operator-= (const value_type& val) {for(std::size_t i = 0; i < N; ++i) m_data[i] -= val;return *this;}
		MathVector& operator*= (const value_type& val) {for(std::size_t i = 0; i < N; ++i) m_data[i] *= val;return *this;}
		MathVector& operator/= (const value_type& val) {for(std::size_t i = 0; i < N; ++i) m_data[i] /= val;return *this;}

		// negation
		MathVector operator- () const { MathVector<N, T> v; for (std::size_t i = 0; i < N; ++i) v.set_coord(i, -m_data[i]); return v; }

		// scalar product
		value_type operator* (const MathVector& v) const {value_type res = 0.0; for(std::size_t i = 0; i < N; ++i) res += m_data[i] * v.coord(i);return res;}

		inline std::size_t size() const {return N;}

		inline value_type& coord(size_t index)				{return m_data[index];}
		inline value_type coord(size_t index) const			{return m_data[index];}
		
		inline value_type& operator[](size_t index)				{return m_data[index];}
		inline const value_type& operator[](size_t index) const	{return m_data[index];}

		inline void set_coord(std::size_t index, value_type v)	{m_data[index] = v;}

	protected:
		value_type	m_data[N];

	protected:
		inline void assign(const MathVector<N>& v) {for(std::size_t i = 0; i < N; ++i) m_data[i] = v.coord(i);}

};

/** THIS IS A TEMPORARY STUB. Some discretization routines use MathVector<0>, requiring
 * a single coordinate.
 * A mathematical Vector with 1 entry and static storage
 */
template <typename T>
class MathVector<0, T>
{
	public:
		typedef std::size_t size_type;
		typedef T value_type;
		static const std::size_t Size = 1;

	public:
		MathVector()	{}
		MathVector(value_type x)
		{
			m_data[0] = x;
		}
		MathVector(const MathVector<0, T>& v)	{assign(v);}

		static inline MathVector from(const MathVector<0, T>& v)	{return v;}
		static inline MathVector from(const MathVector<1, T>& v)	{return MathVector();}
		static inline MathVector from(const MathVector<2, T>& v)	{return MathVector();}
		static inline MathVector from(const MathVector<3, T>& v)	{return MathVector();}
		static inline MathVector from(const MathVector<4, T>& v)	{return MathVector();}
		template <std::size_t fromN>
		static inline MathVector from(const MathVector<fromN, T>& v)
		{
			return MathVectorFrom<0, fromN, T>(v);
		}

		// operations with other vectors
		MathVector& operator=  (const MathVector& v) {assign(v); return *this;}
		MathVector& operator+= (const MathVector& v) {m_data[0] += v.x(); return *this;}
		MathVector& operator-= (const MathVector& v) {m_data[0] -= v.x(); return *this;}

		// operations with scalar
		MathVector& operator=  (const value_type& val) {m_data[0] =  val;return *this;}
		MathVector& operator+= (const value_type& val) {m_data[0] += val;return *this;}
		MathVector& operator-= (const value_type& val) {m_data[0] -= val;return *this;}
		MathVector& operator*= (const value_type& val) {m_data[0] *= val;return *this;}
		MathVector& operator/= (const value_type& val) {m_data[0] /= val;return *this;}

		// negation
		MathVector& operator- () { return MathVector<0, T>(-m_data[0]); }

		// scalar product
		value_type operator* (const MathVector& v) const {return m_data[0] * v.x();}

		inline std::size_t size() const								{return 1;}

		inline value_type& coord(std::size_t index)					{return m_data[0];}
		inline value_type coord(std::size_t index) const			{return m_data[0];}

		inline value_type& operator[](std::size_t index)				{return m_data[0];}
		inline const value_type& operator[](std::size_t index) const	{return m_data[0];}

		inline void set_coord(std::size_t index, value_type v)	{m_data[0] = v;}

		inline value_type& x()						{return m_data[0];}
		inline const value_type& x() const			{return m_data[0];}

		value_type m_data[1];
	protected:
		inline void assign(const MathVector<0, T>& v)	{m_data[0] = v.m_data[0];}

};

/**
 * A mathematical Vector with 1 entry and static storage
 */
template <typename T>
class MathVector<1, T>
{
	public:
		typedef std::size_t size_type;
		typedef T value_type;
		static const std::size_t Size = 1;

	public:
		MathVector() 	{m_data[0] = 0.0;}
		MathVector(value_type x) { m_data[0] = x; }
		MathVector(const MathVector<1, T>& v)	{assign(v);}

		static inline MathVector from(const MathVector<0, T>& v)	{return MathVector(0);}
		static inline MathVector from(const MathVector<1, T>& v)	{return v;}
		static inline MathVector from(const MathVector<2, T>& v)	{return MathVector(v[0]);}
		static inline MathVector from(const MathVector<3, T>& v)	{return MathVector(v[0]);}
		static inline MathVector from(const MathVector<4, T>& v)	{return MathVector(v[0]);}
		template <std::size_t fromN>
		static inline MathVector from(const MathVector<fromN, T>& v)
		{
			return MathVectorFrom<1, fromN, T>(v);
		}

		// operations with other vectors
		MathVector& operator=  (const MathVector& v) {assign(v); return *this;}
		MathVector& operator+= (const MathVector& v) {m_data[0] += v.x(); return *this;}
		MathVector& operator-= (const MathVector& v) {m_data[0] -= v.x(); return *this;}

		// operations with scalar
		MathVector& operator=  (const value_type& val) {m_data[0] =  val;return *this;}
		MathVector& operator+= (const value_type& val) {m_data[0] += val;return *this;}
		MathVector& operator-= (const value_type& val) {m_data[0] -= val;return *this;}
		MathVector& operator*= (const value_type& val) {m_data[0] *= val;return *this;}
		MathVector& operator/= (const value_type& val) {m_data[0] /= val;return *this;}

		// negation
		MathVector operator- () const { return MathVector<1, T>(-m_data[0]); }

		// scalar product
		value_type operator* (const MathVector& v) const {return m_data[0] * v.x();}

		inline std::size_t size() const								{return 1;}

		inline value_type& coord(std::size_t index)					{return m_data[0];}
		inline value_type coord(std::size_t index) const			{return m_data[0];}

		inline value_type& operator[](std::size_t index)				{return m_data[0];}
		inline const value_type& operator[](std::size_t index) const	{return m_data[0];}

		inline void set_coord(std::size_t index, value_type v)	{m_data[0] = v;}

		inline value_type& x()						{return m_data[0];}
		inline const value_type& x() const			{return m_data[0];}

		value_type m_data[1];
	protected:
		inline void assign(const MathVector<1, T>& v)	{m_data[0] = v.m_data[0];}

};

/**
 * A mathematical Vector with 2 entries and static storage
 */
template <typename T>
class MathVector<2, T>
{
	public:
		typedef std::size_t size_type;
		typedef T value_type;
		static const std::size_t Size = 2;

	public:
		MathVector()	{m_data[0] = m_data[1] = 0.0;}
		MathVector(const value_type& val) {m_data[0] = m_data[1] = val;}
		MathVector(value_type x, value_type y)
		{
			m_data[0] = x;
			m_data[1] = y;
		}
		MathVector(const MathVector<2,T>& v)	{assign(v);}

		static inline MathVector from(const MathVector<0, T>& v)	{return MathVector(0, 0);}
		static inline MathVector from(const MathVector<1, T>& v)	{return MathVector(v[0], 0);}
		static inline MathVector from(const MathVector<2, T>& v)	{return v;}
		static inline MathVector from(const MathVector<3, T>& v)	{return MathVector(v[0], v[1]);}
		static inline MathVector from(const MathVector<4, T>& v)	{return MathVector(v[0], v[1]);}
		template <std::size_t fromN>
		static inline MathVector from(const MathVector<fromN, T>& v)
		{
			return MathVectorFrom<2, fromN, T>(v);
		}

		// operations with other vectors
		MathVector& operator=  (const MathVector& v) {assign(v); return *this;}
		MathVector& operator+= (const MathVector& v) {for(std::size_t i = 0; i < 2; ++i) m_data[i] += v.coord(i);return *this;}
		MathVector& operator-= (const MathVector& v) {for(std::size_t i = 0; i < 2; ++i) m_data[i] -= v.coord(i);return *this;}

		// operations with scalar
		MathVector& operator=  (const value_type& val) {for(std::size_t i = 0; i < 2; ++i) m_data[i] =  val;return *this;}
		MathVector& operator+= (const value_type& val) {for(std::size_t i = 0; i < 2; ++i) m_data[i] += val;return *this;}
		MathVector& operator-= (const value_type& val) {for(std::size_t i = 0; i < 2; ++i) m_data[i] -= val;return *this;}
		MathVector& operator*= (const value_type& val) {for(std::size_t i = 0; i < 2; ++i) m_data[i] *= val;return *this;}
		MathVector& operator/= (const value_type& val) {for(std::size_t i = 0; i < 2; ++i) m_data[i] /= val;return *this;}

		// negation
		MathVector operator- () const { return MathVector<2, T>(-m_data[0], -m_data[1]); }

		// scalar product
		value_type operator* (const MathVector& v) const {value_type res = 0.0; for(std::size_t i = 0; i < 2; ++i) res += m_data[i] * v.coord(i);return res;}

		inline std::size_t size() const								{return 2;}

		inline value_type& coord(std::size_t index)					{return m_data[index];}
		inline value_type coord(std::size_t index) const			{return m_data[index];}

		inline value_type& operator[](std::size_t index)				{return m_data[index];}
		inline const value_type& operator[](std::size_t index) const	{return m_data[index];}

		inline void set_coord(std::size_t index, value_type v)	{m_data[index] = v;}

		inline value_type& x()						{return m_data[0];}
		inline const value_type& x() const			{return m_data[0];}

		inline value_type& y()						{return m_data[1];}
		inline const value_type& y() const			{return m_data[1];}

		value_type m_data[2];
	protected:
		inline void assign(const MathVector<2,T>& v)	{m_data[0] = v.m_data[0];
														 m_data[1] = v.m_data[1];}
};

/**
 * A mathematical Vector with 3 entries and static storage
 */
template <typename T>
class MathVector<3, T>
{
	public:
		typedef std::size_t size_type;
		typedef T value_type;
		static const std::size_t Size = 3;

	public:
		MathVector()	{m_data[0] = m_data[1] = m_data[2] = 0.0;}
		MathVector(const value_type& val) {m_data[0] = m_data[1] = m_data[2] = val;}
		MathVector(value_type x, value_type y, value_type z)
		{
			m_data[0] = x;
			m_data[1] = y;
			m_data[2] = z;
		}
		MathVector(const MathVector<3,T>& v)	{assign(v);}

		static inline MathVector from(const MathVector<0, T>& v)	{return MathVector(0, 0, 0);}
		static inline MathVector from(const MathVector<1, T>& v)	{return MathVector(v[0], 0, 0);}
		static inline MathVector from(const MathVector<2, T>& v)	{return MathVector(v[0], v[1], 0);}
		static inline MathVector from(const MathVector<3, T>& v)	{return v;}
		static inline MathVector from(const MathVector<4, T>& v)	{return MathVector(v[0], v[1], v[2]);}
		template <std::size_t fromN>
		static inline MathVector from(const MathVector<fromN, T>& v)
		{
			return MathVectorFrom<3, fromN, T>(v);
		}

		// operations with other vectors
		MathVector& operator=  (const MathVector& v) {assign(v); return *this;}
		MathVector& operator+= (const MathVector& v) {for(std::size_t i = 0; i < 3; ++i) m_data[i] += v.coord(i);return *this;}
		MathVector& operator-= (const MathVector& v) {for(std::size_t i = 0; i < 3; ++i) m_data[i] -= v.coord(i);return *this;}

		// operations with scalar
		MathVector& operator=  (const value_type& val) {for(std::size_t i = 0; i < 3; ++i) m_data[i] =  val;return *this;}
		MathVector& operator+= (const value_type& val) {for(std::size_t i = 0; i < 3; ++i) m_data[i] += val;return *this;}
		MathVector& operator-= (const value_type& val) {for(std::size_t i = 0; i < 3; ++i) m_data[i] -= val;return *this;}
		MathVector& operator*= (const value_type& val) {for(std::size_t i = 0; i < 3; ++i) m_data[i] *= val;return *this;}
		MathVector& operator/= (const value_type& val) {for(std::size_t i = 0; i < 3; ++i) m_data[i] /= val;return *this;}

		// negation
		MathVector operator- () const { return MathVector<3, T>(-m_data[0], -m_data[1], -m_data[2]); }

		// scalar product
		value_type operator* (const MathVector& v) const {value_type res = 0.0; for(std::size_t i = 0; i < 3; ++i) res += m_data[i] * v.coord(i);return res;}

		inline std::size_t size() const										{return 3;}

		inline value_type& coord(std::size_t index)					{return m_data[index];}
		inline value_type coord(std::size_t index) const			{return m_data[index];}

		inline value_type& operator[](std::size_t index)				{return m_data[index];}
		inline const value_type& operator[](std::size_t index) const	{return m_data[index];}

		inline void set_coord(std::size_t index, value_type v)	{m_data[index] = v;}

		inline value_type& x()						{return m_data[0];}
		inline const value_type& x() const			{return m_data[0];}

		inline value_type& y()						{return m_data[1];}
		inline const value_type& y() const			{return m_data[1];}

		inline value_type& z()						{return m_data[2];}
		inline const value_type& z() const			{return m_data[2];}

		value_type m_data[3];
	protected:
		inline void assign(const MathVector<3,T>& v)	{m_data[0] = v.m_data[0];
												 		 m_data[1] = v.m_data[1];
												 		 m_data[2] = v.m_data[2];}

};

/**
 * A mathematical Vector with 4 entries and static storage
 */
template <typename T>
class MathVector<4, T>
{
	public:
		typedef std::size_t size_type;
		typedef T value_type;
		static const std::size_t Size = 4;

	public:
		MathVector()	{m_data[0] = m_data[1] = m_data[2] = m_data[3] = 0.0;}
		MathVector(const value_type& val) {m_data[0] = m_data[1] = m_data[2] = m_data[3] =val;}
		MathVector(value_type x, value_type y, value_type z, value_type w)
		{
			m_data[0] = x;
			m_data[1] = y;
			m_data[2] = z;
			m_data[3] = w;
		}
		MathVector(const MathVector<4,T>& v)	{assign(v);}

		static inline MathVector from(const MathVector<0, T>& v)	{return MathVector(0, 0, 0, 0);}
		static inline MathVector from(const MathVector<1, T>& v)	{return MathVector(v[0], 0, 0, 0);}
		static inline MathVector from(const MathVector<2, T>& v)	{return MathVector(v[0], v[1], 0, 0);}
		static inline MathVector from(const MathVector<3, T>& v)	{return MathVector(v[0], v[1], v[2], 0);}
		static inline MathVector from(const MathVector<4, T>& v)	{return v;}
		template <std::size_t fromN>
		static inline MathVector from(const MathVector<fromN, T>& v)
		{
			return MathVectorFrom<4, fromN, T>(v);
		}

		// operations with other vectors
		MathVector& operator=  (const MathVector& v) {assign(v); return *this;}
		MathVector& operator+= (const MathVector& v) {for(std::size_t i = 0; i < 4; ++i) m_data[i] += v.coord(i);return *this;}
		MathVector& operator-= (const MathVector& v) {for(std::size_t i = 0; i < 4; ++i) m_data[i] -= v.coord(i);return *this;}

		// operations with scalar
		MathVector& operator=  (const value_type& val) {for(std::size_t i = 0; i < 4; ++i) m_data[i] =  val;return *this;}
		MathVector& operator+= (const value_type& val) {for(std::size_t i = 0; i < 4; ++i) m_data[i] += val;return *this;}
		MathVector& operator-= (const value_type& val) {for(std::size_t i = 0; i < 4; ++i) m_data[i] -= val;return *this;}
		MathVector& operator*= (const value_type& val) {for(std::size_t i = 0; i < 4; ++i) m_data[i] *= val;return *this;}
		MathVector& operator/= (const value_type& val) {for(std::size_t i = 0; i < 4; ++i) m_data[i] /= val;return *this;}

		// negation
		MathVector operator- () const { return MathVector<4, T>(-m_data[0], -m_data[1], -m_data[2], -m_data[3]); }

		// scalar product
		value_type operator* (const MathVector& v) const {value_type res = 0.0; for(std::size_t i = 0; i < 4; ++i) res += m_data[i] * v.coord(i);return res;}

		inline std::size_t size() const									{return 4;}

		inline value_type& coord(std::size_t index)					{return m_data[index];}
		inline value_type coord(std::size_t index) const			{return m_data[index];}

		inline value_type& operator[](std::size_t index)				{return m_data[index];}
		inline const value_type& operator[](std::size_t index) const	{return m_data[index];}

		inline void set_coord(std::size_t index, value_type v)	{m_data[index] = v;}

		inline value_type& x()						{return m_data[0];}
		inline const value_type& x() const			{return m_data[0];}

		inline value_type& y()						{return m_data[1];}
		inline const value_type& y() const			{return m_data[1];}

		inline value_type& z()						{return m_data[2];}
		inline const value_type& z() const			{return m_data[2];}

		inline value_type& w()						{return m_data[3];}
		inline const value_type& w() const			{return m_data[3];}

		value_type m_data[4];
	protected:
		inline void assign(const MathVector<4,T>& v)	{m_data[0] = v.m_data[0];
														 m_data[1] = v.m_data[1];
														 m_data[2] = v.m_data[2];
														 m_data[3] = v.m_data[3];}

};

template <std::size_t N, typename T>
bool operator== (const MathVector<N,T>& v, const MathVector<N,T>& w)
{
	for(std::size_t i = 0; i < N; ++i)
	{
		if(v[i] != w[i]) return false;
	}
	return true;
}

template <typename T>
bool operator== (const MathVector<0,T>& v, const MathVector<0,T>& w)
{
	return true;
}

//	NOTE: this implementation determines the state of the relation '<'
//	by considering the FIRST vector-entry, which is not equal
template <std::size_t N, typename T>
bool operator< (const MathVector<N,T>& v, const MathVector<N,T>& w)
{
	for(std::size_t i = 0; i < N; ++i)
	{
		if(v[i] < w[i] - SMALL) 		return true;
		else if(v[i] > w[i] + SMALL)	return false;
	}
	return false;
}

template <std::size_t N, typename T>
bool operator!= (const MathVector<N,T>& v, const MathVector<N,T>& w)
{
	return !(v == w);
}

template <std::size_t N, typename T>
std::ostream& operator<< (std::ostream& outStream, const ug::MathVector<N,T>& v)
{
	for(std::size_t i = 0; i < N; ++i)
		outStream << "[" << i << "]: " << v.coord(i) << std::endl;
	return outStream;
}

template <typename T>
std::ostream& operator<< (std::ostream& outStream, const ug::MathVector<0,T>& v)
{
	outStream << "(empty)";
	return outStream;
}

template <typename T>
std::ostream& operator<< (std::ostream& outStream, const ug::MathVector<1,T>& v)
{
	outStream << "(" << v[0] << ")";
	return outStream;
}
template <typename T>
std::ostream& operator<< (std::ostream& outStream, const ug::MathVector<2,T>& v)
{
	outStream << "("<<v[0]<<", "<<v[1]<<")";
	return outStream;
}
template <typename T>
std::ostream& operator<< (std::ostream& outStream, const ug::MathVector<3,T>& v)
{
	outStream << "("<<v[0]<<", "<<v[1]<<", "<<v[2]<<")";
	return outStream;
}
template <typename T>
std::ostream& operator<< (std::ostream& outStream, const ug::MathVector<4,T>& v)
{
	outStream << "("<<v[0]<<", "<<v[1]<<", "<<v[2]<<", "<<v[3]<<")";
	return outStream;
}

UG_API std::ostream& operator<< (std::ostream& outStream, const ug::MathVector<1>& v);
UG_API std::ostream& operator<< (std::ostream& outStream, const ug::MathVector<2>& v);
UG_API std::ostream& operator<< (std::ostream& outStream, const ug::MathVector<3>& v);
UG_API std::ostream& operator<< (std::ostream& outStream, const ug::MathVector<4>& v);

UG_API std::ostream& write_plain_txt (std::ostream& outStream, const ug::MathVector<1>& v);
UG_API std::ostream& write_plain_txt (std::ostream& outStream, const ug::MathVector<2>& v);
UG_API std::ostream& write_plain_txt (std::ostream& outStream, const ug::MathVector<3>& v);
UG_API std::ostream& write_plain_txt (std::ostream& outStream, const ug::MathVector<4>& v);

UG_API std::istream& read_plain_txt (std::istream& inStream, ug::MathVector<1>& v);
UG_API std::istream& read_plain_txt (std::istream& inStream, ug::MathVector<2>& v);
UG_API std::istream& read_plain_txt (std::istream& inStream, ug::MathVector<3>& v);
UG_API std::istream& read_plain_txt (std::istream& inStream, ug::MathVector<4>& v);

template <class TStream, std::size_t N, class T>
void Serialize(TStream& out, const MathVector<N, T>& val)
{
	out.write((char*)val.m_data, sizeof(T) * N);
}

template <class TStream, std::size_t N, class T>
void Deserialize(TStream& out, MathVector<N, T>& valOut)
{
	out.read((char*)valOut.m_data, sizeof(T) * N);
}

// end group vectors
/// \}

}//	end of namespace


#endif /* __H__COMMON__MATH_MathVector__ */
