/*
 * Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__PCL_reduce_traits
#define __H__PCL_reduce_traits


#include <algorithm>

#include "common/error.h"

#include "pcl_methods.h"

namespace pcl {
	
///	methods defined in those traits are used by ComPol_AttachmentReduce
/**	A default implementation is provided which works for integer types */
template <typename TValue>
struct reduce_traits
{
	using value_t = TValue;
	static value_t min(value_t v1, value_t v2) {return std::min(v1, v2);}
	static value_t max(value_t v1, value_t v2) {return std::max(v1, v2);}
	static value_t sum(value_t v1, value_t v2) {return v1 + v2;}
	static value_t prod(value_t v1, value_t v2) {return v1 * v2;}
	static value_t land(value_t v1, value_t v2) {return v1 && v2;}
	static value_t band(value_t v1, value_t v2) {return v1 & v2;}
	static value_t lor(value_t v1, value_t v2) {return v1 || v2;}
	static value_t bor(value_t v1, value_t v2) {return v1 | v2;}
};

/**	Specialization for float. No band and bor operations are supported. */
template <>
struct reduce_traits<float>
{
	using value_t = float;
	static value_t min(value_t v1, value_t v2) {return std::min(v1, v2);}
	static value_t max(value_t v1, value_t v2) {return std::max(v1, v2);}
	static value_t sum(value_t v1, value_t v2) {return v1 + v2;}
	static value_t prod(value_t v1, value_t v2) {return v1 * v2;}
	static value_t land(value_t v1, value_t v2) {return v1 && v2;}
	static value_t band(value_t v1, value_t v2) {UG_THROW("floats do not support a binary and operation.");}
	static value_t lor(value_t v1, value_t v2) {return v1 || v2;}
	static value_t bor(value_t v1, value_t v2) {UG_THROW("floats do not support a binary or operation.");}
};

/**	Specialization for double. No band and bor operations are supported. */
template <>
struct reduce_traits<double>
{
	using value_t = double;
	static value_t min(value_t v1, value_t v2) {return std::min(v1, v2);}
	static value_t max(value_t v1, value_t v2) {return std::max(v1, v2);}
	static value_t sum(value_t v1, value_t v2) {return v1 + v2;}
	static value_t prod(value_t v1, value_t v2) {return v1 * v2;}
	static value_t land(value_t v1, value_t v2) {return v1 && v2;}
	static value_t band(value_t v1, value_t v2) {UG_THROW("doubles do not support a binary and operation.");}
	static value_t lor(value_t v1, value_t v2) {return v1 || v2;}
	static value_t bor(value_t v1, value_t v2) {UG_THROW("doubles do not support a binary or operation.");}
};


template <typename T>
class Reducer {
public:
	Reducer (ReduceOperation rop)
	{
		if(rop == MPI_MIN) m_op = &reduce_traits<T>::min;
		else if(rop == MPI_MAX) m_op = &reduce_traits<T>::max;
		else if(rop == MPI_SUM) m_op = &reduce_traits<T>::sum;
		else if(rop == MPI_PROD) m_op = &reduce_traits<T>::prod;
		else if(rop == MPI_LAND) m_op = &reduce_traits<T>::land;
		else if(rop == MPI_BAND) m_op = &reduce_traits<T>::band;
		else if(rop == MPI_LOR) m_op = &reduce_traits<T>::lor;
		else if(rop == MPI_BOR) m_op = &reduce_traits<T>::bor;
		else {UG_THROW ("Unsupported reduce operation: " << rop)};
	}

	
	T operator () (T v1, T v2)	{return m_op (v1, v2);}

	T reduce (T v1, T v2)	{return m_op (v1, v2);}

private:
	T (*m_op) (T, T);
};



}//	end of namespace

#endif	//__H__PCL_reduce_traits
