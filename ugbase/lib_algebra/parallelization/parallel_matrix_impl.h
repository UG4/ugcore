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

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX_IMPL__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX_IMPL__

#include "parallel_matrix.h"

namespace ug
{

template <typename TMatrix>
typename ParallelMatrix<TMatrix>::this_type&
ParallelMatrix<TMatrix>::operator =(const typename ParallelMatrix<TMatrix>::this_type &M)
{
//	forward to sequential matrices
	TMatrix::operator= (*dynamic_cast<const TMatrix*>(&M));

//	copy storage type and layouts
	this->set_storage_type(M.get_storage_mask());
	this->set_layouts(M.layouts());

//	we're done
	return *this;
}

template <typename TMatrix>
bool
ParallelMatrix<TMatrix>::
change_storage_type(ParallelStorageType type)
{
	// can only change if current state is defined
	if(has_storage_type(PST_UNDEFINED))
	{
		UG_LOG("Current Storage Type is undefined. "
				"Cannot change storage type.\n");
		return false;
	}

	// if already in that type
	if(has_storage_type(type))
		return true;

//	todo: Implement more types.
	UG_THROW("ParallelMatrix::change_storage_type:"
			" Currently no storage conversion supported.");
}

// calculate res = A x
template <typename TMatrix>
template<typename TPVector>
bool
ParallelMatrix<TMatrix>::
apply(TPVector &res, const TPVector &x) const
{
	PROFILE_FUNC_GROUP("algebra");
//	check types combinations
	int type = -1;
	if(has_storage_type(PST_ADDITIVE)
			&& x.has_storage_type(PST_CONSISTENT)) type = 0;
	if(has_storage_type(PST_CONSISTENT)
			&& x.has_storage_type(PST_ADDITIVE)) type = 1;
	if(has_storage_type(PST_CONSISTENT)
			&& x.has_storage_type(PST_CONSISTENT)) type = 2;

//	if no admissible type is found, return error
	if(type == -1)
	{
		UG_THROW("ParallelMatrix::apply (b = A*x): "
				"Wrong storage type of Matrix/Vector: Possibilities are:\n"
				"    - A is PST_ADDITIVE and x is PST_CONSISTENT\n"
				"    - A is PST_CONSISTENT and x is PST_ADDITIVE\n"
				"    (storage type of A = " << get_storage_type() << ", x = " << x.get_storage_type() << ")");
	}

//	apply on single process vector
	TMatrix::axpy(res, 0.0, res, 1.0, x);

//	set outgoing vector to additive storage
	switch(type)
	{
		case 0: res.set_storage_type(PST_ADDITIVE); break;
		case 1: res.set_storage_type(PST_ADDITIVE); break;
		case 2: res.set_storage_type(PST_CONSISTENT); break;
	}

//	we're done.
	return true;
}

// calculate res = A.T x
template <typename TMatrix>
template<typename TPVector>
bool
ParallelMatrix<TMatrix>::
apply_transposed(TPVector &res, const TPVector &x) const
{
	PROFILE_FUNC_GROUP("algebra");
//	check types combinations
	int type = -1;
	if(has_storage_type(PST_ADDITIVE)
			&& x.has_storage_type(PST_CONSISTENT)) type = 0;
	if(has_storage_type(PST_CONSISTENT)
			&& x.has_storage_type(PST_ADDITIVE)) type = 1;
	if(has_storage_type(PST_CONSISTENT)
			&& x.has_storage_type(PST_CONSISTENT)) type = 2;

//	if no admissible type is found, return error
	if(type == -1)
	{
		UG_THROW("ParallelMatrix::apply_transposed (b = A^T*x): "
				"Wrong storage type of Matrix/Vector: Possibilities are:\n"
				"    - A is PST_ADDITIVE and x is PST_CONSISTENT\n"
				"    - A is PST_CONSISTENT and x is PST_ADDITIVE\n"
				"    (storage type of A = " << get_storage_type() << ", x = " << x.get_storage_type() << ")");
	}

//	apply on single process vector
	TMatrix::axpy_transposed(res, 0.0, res, 1.0, x);

//	set outgoing vector to additive storage
	switch(type)
	{
		case 0: res.set_storage_type(PST_ADDITIVE); break;
		case 1: res.set_storage_type(PST_ADDITIVE); break;
		case 2: res.set_storage_type(PST_CONSISTENT); break;
	}

//	we're done.
	return true;
}

// calculate res -= A x
template <typename TMatrix>
template<typename TPVector>
bool
ParallelMatrix<TMatrix>::
matmul_minus(TPVector &res, const TPVector &x) const
{
	PROFILE_FUNC_GROUP("algebra");
//	check types combinations
	int type = -1;
	if(this->has_storage_type(PST_ADDITIVE)
			&& x.has_storage_type(PST_CONSISTENT)
			&& res.has_storage_type(PST_ADDITIVE)) type = 0;

//	if no admissible type is found, return error
	if(type == -1)
	{
		UG_THROW("ParallelMatrix::matmul_minus (b -= A*x):"
				" Wrong storage type of Matrix/Vector: Possibilities are:\n"
				"    - A is PST_ADDITIVE and x is PST_CONSISTENT and b is PST_ADDITIVE\n"
				"    (storage type of A = " << this->get_storage_type() << ", x = " << x.get_storage_type() << ", b = " << res.get_storage_type() << ")");
	}

//	apply on single process vector
	TMatrix::axpy(res, 1.0, res, -1.0, x);

//	set outgoing vector to additive storage
//	(it could have been PST_UNIQUE before)
	switch(type)
	{
		case 0: res.set_storage_type(PST_ADDITIVE); break;
	}

//	we're done.
	return true;
}


template<typename matrix_type, typename vector_type>
ParallelStorageType GetMultType(const ParallelMatrix<matrix_type> &A1, const ParallelVector<vector_type> &x)
{
	ParallelStorageType type = PST_UNDEFINED;
	if(A1.has_storage_type(PST_ADDITIVE)
			&& x.has_storage_type(PST_CONSISTENT)) type = PST_ADDITIVE; // 0
	if(A1.has_storage_type(PST_CONSISTENT)
			&& x.has_storage_type(PST_ADDITIVE)) type = PST_ADDITIVE; // 1
	if(A1.has_storage_type(PST_CONSISTENT)
			&& x.has_storage_type(PST_CONSISTENT)) type = PST_CONSISTENT; // 2

	//	if no admissible type is found, return error
	if(type == PST_UNDEFINED)
	{
		UG_THROW("ParallelMatrix::apply (b = A*x): "
				"Wrong storage type of Matrix/Vector: Possibilities are:\n"
				"    - A is PST_ADDITIVE and x is PST_CONSISTENT\n"
				"    - A is PST_CONSISTENT and x is PST_ADDITIVE\n"
				"    (storage type of A = " << A1.get_storage_type() << ", x = " << x.get_storage_type() << ")");
	}
	return type;
}

template<typename matrix_type, typename vector_type>
inline bool MatMultDirect(ParallelVector<vector_type> &dest,
		const number &beta1, const ParallelMatrix<matrix_type> &A1, const ParallelVector<vector_type> &w1)
{
	PROFILE_FUNC_GROUP("algebra");
	//	check types combinations
	ParallelStorageType type = GetMultType(A1, w1);
	if(type == PST_UNDEFINED) return false;

	MatMult(dynamic_cast<vector_type&>(dest), beta1, dynamic_cast<const matrix_type&>(A1), dynamic_cast<const vector_type&>(w1));

	//	set outgoing vector storage type
	dest.set_storage_type(type);

	//	we're done.
	return true;
}

template<typename matrix_type, typename vector_type>
inline bool MatMultAddDirect(ParallelVector<vector_type> &dest,
		const number &alpha1, const ParallelVector<vector_type> &v1,
		const number &beta1, const ParallelMatrix<matrix_type> &A1, const ParallelVector<vector_type> &w1)
{
	PROFILE_FUNC_GROUP("algebra");
	//	check types combinations
	ParallelStorageType type = GetMultType(A1, w1);
	if(type == PST_UNDEFINED) return false;

	if(!v1.has_storage_type(type))
	{
		UG_THROW("MatMultAdd(dest, alpha1, v1, beta1, A1, w1): Storage type of A1*w1 = " << type << ", storage type of v1 = " << v1.get_storage_type() << ".");
	}

	MatMultAdd(dynamic_cast<vector_type&>(dest), alpha1, dynamic_cast<const vector_type&>(v1),
			beta1, dynamic_cast<const matrix_type&>(A1), dynamic_cast<const vector_type&>(w1));

	//	set outgoing vector storage type
	dest.set_storage_type(type);

	//	we're done.
	return true;
}

template<typename matrix_type, typename vector_type>
inline bool MatMultAddDirect(ParallelVector<vector_type> &dest,
		const number &alpha1, const ParallelVector<vector_type> &v1,
		const number &alpha2, const ParallelVector<vector_type> &v2,
		const number &beta1, const ParallelMatrix<matrix_type> &A1, const ParallelVector<vector_type> &w1)
{
	PROFILE_FUNC_GROUP("algebra");
	//	check types combinations
	ParallelStorageType type = GetMultType(A1, w1);
	if(type == PST_UNDEFINED) return false;

	if(!v1.has_storage_type(type) || !v2.has_storage_type(type))
	{
		UG_THROW("MatMultAdd(dest, alpha1, v1, alpha2, v2, beta1, A1, w1):"
				" Storage type of A1*w1 doesnt match storage type of v1 or v2.\n"
				"Storage type of A1*w1 = " << type << " v1 = " << v1.get_storage_type() << ", v2 = " << v2.get_storage_type() << ".");
	}

	MatMultAdd(dynamic_cast<vector_type&>(dest), alpha1, dynamic_cast<const vector_type&>(v1), alpha2, dynamic_cast<const vector_type&>(v2),
			beta1, dynamic_cast<const matrix_type&>(A1), dynamic_cast<const vector_type&>(w1));

	//	set outgoing vector storage type
	dest.set_storage_type(type);

	//	we're done.
	return true;
}

template<typename matrix_type, typename vector_type>
inline void MatMultTransposedDirect(ParallelVector<vector_type> &dest,
		const number &beta1, const ParallelMatrix<matrix_type> &A1, const ParallelVector<vector_type> &w1)
{
	PROFILE_FUNC_GROUP("algebra");
	//	check types combinations
	ParallelStorageType type = GetMultType(A1, w1);
	if(type == PST_UNDEFINED) UG_THROW("Wrong storage type");

	MatMultTransposed(dynamic_cast<vector_type&>(dest), beta1, dynamic_cast<const matrix_type&>(A1), dynamic_cast<const vector_type&>(w1));

	//	set outgoing vector storage type
	dest.set_storage_type(type);
}


} // end namespace ug

#endif