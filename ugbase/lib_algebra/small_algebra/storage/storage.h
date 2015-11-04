/*
 *  storage.h
 *
 *  Created by Martin Rupp on 21.07.10.
 *  Copyright 2010 . All rights reserved.
 *
 */

#ifndef __H__UG__COMMON__STORAGE_H__
#define __H__UG__COMMON__STORAGE_H__


namespace ug{

/**
 * Use eMatrixOrdering to indicate if in a two-dimensional array
 * RowMajor: each row is stored consecutively in memory, or
 * ColMajor: each col is stored consecutively in memory.
 * \sa http://en.wikipedia.org/wiki/Row-major_order .
  */
enum eMatrixOrdering
{
	RowMajor = 0,
	ColMajor = 1
};

/**
 * use these traits to distinguish between static and variable types
 * like FixedArray1, VariableArray1 and stl::vector.
 */
template<typename TStorage>
struct storage_traits1
{
	enum {is_static = false};
	enum {static_size = 0};
};

struct static_type{};
struct variable_type{};



}

#include "fixed_array.h"
#include "variable_array.h"

#endif // __H__UG__COMMON__MATRIX_FLAGS_H__
