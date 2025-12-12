/*
 * Copyright (c) 2010-2011:  G-CSC, Goethe University Frankfurt
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


#ifndef __H__UG__SMALL_ALGEBRA__BLOCKS__
#define __H__UG__SMALL_ALGEBRA__BLOCKS__


namespace ug {
	
inline double dabs(double a) { return a > 0 ? a : -a; }


template <typename t> struct block_traits;
template<typename value_type, typename vec_type> struct block_multiply_traits;


//////////////////////////////////////////////////////

template<typename TYPE>
inline double BlockNorm2(const TYPE &v)
{
	return v.norm2();
}

template<typename TYPE>
inline double BlockNorm(const TYPE &v)
{
	return sqrt(BlockNorm2(v));
}


//////////////////////////////////////////////////////

// get/set vector
template<typename T> inline double &BlockRef(T &vec, size_t i)
{
	return vec[i];
}

template<typename T> inline const double &BlockRef(const T &vec, size_t i)
{
	return vec[i];
}

// get/set matrix
template<typename T> inline double &BlockRef(T &mat, size_t i, size_t j)
{
	return mat(i, j);
}

template<typename T> inline const double &BlockRef(const T &mat, size_t i, size_t j)
{
	return mat(i, j);
}


//////////////////////////////////////////////////////

/*template<typename T, size_t n> inline bool BlockSerialize(const DenseVector<FixedArray<T, n> > &t, std::ostream &buff)
{
	for(int i=0; i<n; i++)
		BlockSerialize(v[i], buff);
	return true;
}


template<typename T, size_t n> inline bool BlockDeserialize(const DenseVector<FixedArray<T, n> > &t, std::ostream &buff)
{
	for(int i=0; i<n; i++)
		BlockDeserialize(v[i], buff);
	return true;
}
*/
//////////////////////////////////////////////////////


// algebra stuff to avoid temporary variables 

	
// MATRICES

// todo: replace add_mult etc. with template expressions
// dest = b*vec
template<typename A, typename B, typename C> inline void AssignMult(A &dest, const B &b, const C &vec);
// dest += b*vec
template<typename A, typename B, typename C> inline void AddMult(A &dest, const B &b, const C &vec);
// dest -= b*vec
template<typename A, typename B, typename C> inline void SubMult(A &dest, const B &b, const C &vec);

// VECTORs



//////////////////////////////////////////////////////
//setSize(t, a, b) for matrices
template<typename T>
inline void SetSize(T &t, size_t a, size_t b);

//setSize(t, a) for vectors
template<typename T>
inline void SetSize(T &t, size_t a);

// getSize
template<typename T>
inline size_t GetSize(const T &t);

//getRows
template<typename T>
inline size_t GetRows(const T &t);

// getRows
template<typename T>
inline size_t GetCols(const T &t);



} // namespace ug

#include "double.h"
#include "small_matrix/densevector.h"
#include "small_matrix/densematrix.h"
#include "small_matrix/block_dense.h"
#include "storage/storage.h"

#endif