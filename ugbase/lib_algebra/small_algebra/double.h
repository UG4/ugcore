/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

/*
 *  blockMatrix.h
 *
 *  Header File for general block matrix / double accessing 
 *  i.e. setAt(mat, i, j, d) -> mat(i,j) = d
 *	and setAt(f, 0, 0, d) -> f = d.
 *  This means doubles and block matrices can be accessed
 *	by the same methods.
 */

// todo: also with float / complex<float> / complex<double>

#ifndef __H__UG__SMALL_ALGEBRA__DOUBLE__
#define __H__UG__SMALL_ALGEBRA__DOUBLE__

#include "blocks.h"
#include "common/common.h"

namespace ug{


//////////////////////////////////////////////////////
template<typename T> number BlockNorm(const T &t); // ø todo check if declaration is redundand
template <>
inline number BlockNorm(const number &t)
{
	return t>0 ? t : -t;
}

template<typename T> number BlockNorm2(const T &t); // ø todo check if declaration is redundand
template <>
inline number BlockNorm2(const number &a)
{
	return a*a;
}

template<typename T> number BlockMaxNorm(const T &t); // ø todo check if declaration is redundand
template <>
inline number BlockMaxNorm(const number &a)
{
	return a>0 ? a : -a;
}

//////////////////////////////////////////////////////
// get/set specialization for numbers

template<> inline number &BlockRef(number &m, size_t i)
{
	UG_ASSERT(i == 0, "block is number, doesnt have component (" << i << ").");
	return m;
}
template<> inline const number &BlockRef(const number &m, size_t i)
{
	UG_ASSERT(i == 0, "block is number, doesnt have component (" << i << ").");
	return m;
}

template<> inline number &BlockRef(number &m, size_t i, size_t j)
{
	UG_ASSERT(i == 0 && j == 0, "block is number, doesnt have component (" << i << ", " << j << ").");
	return m;
}
template<> inline const number &BlockRef(const number &m, size_t i, size_t j)
{
	UG_ASSERT(i == 0 && j == 0, "block is number, doesnt have component (" << i << ", " << j << ").");
	return m;
}

//////////////////////////////////////////////////////
// algebra stuff to avoid temporary variables 

inline void AssignMult(number &dest, const number &b, const number &vec)
{
	dest = b*vec;
}
// dest += vec*b
inline void AddMult(number &dest, const number &b, const number &vec)
{
	dest += b*vec;
}


// dest -= vec*b
inline void SubMult(number &dest, const number &b, const number &vec)
{
	dest -= b*vec;
}


//////////////////////////////////////////////////////
//setSize(t, a, b) for numbers
template<>
inline void SetSize(number &d, size_t a)
{
	UG_ASSERT(a == 1, "block is number, cannot change size to " << a << ".");
}

template<>
inline void SetSize(number &d, size_t a, size_t b)
{
	UG_ASSERT(a == 1 && b == 1, "block is number, cannot change size to (" << a << ", " << b << ").");
}

template<>
inline size_t GetSize(const number &t)
{
	return 1;
}

template<>
inline size_t GetRows(const number &t)
{
	return 1;
}

template<>
inline size_t GetCols(const number &t)
{
	return 1;
}
///////////////////////////////////////////////////////////////////

inline bool InverseMatMult(number &dest, const double &beta, const number &mat, const number &vec)
{
	dest = beta*vec/mat;
	return true;
}

///////////////////////////////////////////////////////////////////
// traits: information for numbers


template<>
struct block_traits<number>
{
	using vec_type = number;
	using inverse_type = number;

	enum { is_static = true};
	enum { static_num_rows = 1};
	enum { static_num_cols = 1};
	enum { static_size = 1 };
	enum { depth = 0 };
};

template<> struct block_multiply_traits<number, number>
{
	using ReturnType = number;
};

template<typename T>
struct block_multiply_traits<number, T>
{
	using ReturnType = T;
};

template<typename T>
struct block_multiply_traits<T, number>
{
	using ReturnType = T;
};

inline bool GetInverse(number &inv, const number &m)
{
	inv = 1.0/m;
	return (m != 0.0);
}

inline bool Invert(number &m)
{
	bool b = (m != 0.0);
	m = 1 / m;
	return b;
}

} // namespace ug

#endif



