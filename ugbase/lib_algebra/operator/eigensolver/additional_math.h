/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __UG__ADDITIONAL_MATH_H__
#define __UG__ADDITIONAL_MATH_H__

#include "smart_ptr_vector.h"

#define PINVIT_PROFILE_FUNC() PROFILE_FUNC_GROUP("pinvit algebra")
#define PINVIT_PROFILE_BEGIN(t) PROFILE_BEGIN_GROUP(t, "pinvit algebra")
#define PINVIT_PROFILE_END() PROFILE_END()

namespace ug{

/*template<typename mat_type, typename vec_type, typename densematrix_type>
void MultiEnergyProd(const SparseMatrix<mat_type> &A,
			Vector<vec_type> **x,
			DenseMatrix<densematrix_type> &rA, size_t n)
{
	UG_ASSERT(n == rA.num_rows() && n == rA.num_cols(), "");
	vec_type Ai_xc;
	rA = 0.0;
	for(size_t i=0; i<A.num_rows(); i++)
	{
		for(size_t c=0; c<n; c++)
		{
			// Ai_xc = A[i] * x[c].
			Ai_xc = 0.0;
			A.mat_mult_add_row(i, Ai_xc, 1.0, (*x[c]));
			for(size_t r=0; r<n; r++)
				rA(c, r) += VecDot((*x[r])[i], Ai_xc);
		}
	}
}*/


inline bool absCompare(double a, double b)
{
	return fabs(a) < fabs(b);
}



template<typename vector_type, typename densematrix_type>
void MultiScalProd(vector_type &px,
			DenseMatrix<densematrix_type> &rA, size_t n)
{
	PINVIT_PROFILE_FUNC();
//	UG_ASSERT(0, "");
	UG_ASSERT(n == rA.num_rows() && n == rA.num_cols(), "");
	for(size_t r=0; r<n; r++)
		for(size_t c=r; c<n; c++)
			rA(r, c) = px[c]->dotprod(*px[r]);

	for(size_t r=0; r<n; r++)
		for(size_t c=0; c<r; c++)
			rA(r,c) = rA(c, r);
}

template<typename matrix_type, typename vector_type>
double EnergyProd(vector_type &v1, matrix_type &A, vector_type &v2, vector_type &tmp)
{
	PINVIT_PROFILE_FUNC();

#ifdef UG_PARALLEL
	pcl::ProcessCommunicator pc;
	v2.change_storage_type(PST_CONSISTENT);
#endif
	A.apply(tmp, v2);
	// tmp is additive, v1 is consistent
	double a = v1.dotprod(tmp);
	//UG_LOG("EnergyProd " << a << "\n");

	return a;
}

template<typename matrix_type, typename vector_type>
double EnergyProd(vector_type &v1, matrix_type &A, vector_type &v2)
{
	vector_type t;
	CloneVector(t, v1);
	return EnergyProd(v1, A, v2, t);
}

template<typename matrix_type, typename vector_type>
double EnergyNorm(vector_type &x, matrix_type &A, vector_type &tmp)
{
	return sqrt( EnergyProd(x, A, x, tmp) );
}

template<typename matrix_type, typename vector_type>
double EnergyNorm(vector_type &x, matrix_type &A)
{
	vector_type tmp;
	CloneVector(tmp, x);
	return sqrt( EnergyProd(x, A, x, tmp) );
}


template<typename matrix_type, typename vector_type, typename densematrix_type>
void MultiEnergyProd(matrix_type &A,
			SmartPtrVector<vector_type> &px,
			DenseMatrix<densematrix_type> &rA, size_t n)
{
	PINVIT_PROFILE_FUNC();
#ifdef UG_PARALLEL
	pcl::ProcessCommunicator pc;
#endif
	UG_ASSERT(n == rA.num_rows() && n == rA.num_cols(), "");
	vector_type t;

#if 0
	CloneVector(t, px(0));

	for(size_t r=0; r<n; r++)
	{
//		t.set(0.0);
	#ifdef UG_PARALLEL
		px(r).change_storage_type(PST_CONSISTENT);
	#endif
		A.apply(t, px(r));
		// tmp is additive, v1 is consistent
		//UG_LOG("EnergyProd " << a << "\n");
		for(size_t c=r; c<n; c++)
		{
			double a = px(c).dotprod(t);
			rA(c, r) = rA(r, c) = a; //EnergyProd(px(r), A, px(c), t);
		}
	}

#else
	CloneVector(t, *px[0]);

#ifdef UG_PARALLEL
	for(size_t r=0; r<n; r++)
		px[r]->change_storage_type(PST_CONSISTENT);
#endif

	for(size_t r=0; r<n; r++)
	{
		// todo: why is SparseMatrix<T>::apply not const ?!?
		A.apply(t, *px[r]);
		// t additive

		for(size_t c=r; c<n; c++)
			rA(c, r) = rA(r, c) = px[c]->dotprod(t);
	}

#endif
}


template<typename tvector>
void PrintStorageType(const tvector &v)
{
#ifdef UG_PARALLEL
	if(v.has_storage_type(PST_UNDEFINED))
		UG_LOG("PST_UNDEFINED ");
	if(v.has_storage_type(PST_CONSISTENT))
		UG_LOG("PST_CONSISTENT ");
	if(v.has_storage_type(PST_ADDITIVE))
		UG_LOG("PST_ADDITIVE ");
	if(v.has_storage_type(PST_UNIQUE))
		UG_LOG("PST_UNIQUE ");
#else
	UG_LOG("serial ");
#endif
}


template<typename matrix_type>
void PrintMatrix(const matrix_type &mat, const char *name)
{
	UG_LOG(name << ":\n" << name << " := matrix([\n");
	for(size_t r=0; r<mat.num_rows(); r++)
	{
		UG_LOG("[");
		for(size_t c=0; c<mat.num_cols(); c++)
		{
			UG_LOG(mat(r, c));
			if(c < mat.num_cols()-1) UG_LOG(",\t");
		}
		UG_LOG("]\n");
	}
	UG_LOG("]);\n");

}

template<typename matrix_type>
void PrintMaple(const matrix_type &mat, const char *name)
{
	UG_LOG(name << ":\n" << name << " := matrix([");
	for(size_t r=0; r<mat.num_rows(); r++)
	{
		UG_LOG("[");
		for(size_t c=0; c<mat.num_cols(); c++)
		{
			UG_LOG(mat(r, c));
			if(c < mat.num_cols()-1) UG_LOG(", ");
		}
		UG_LOG("]");
		if(r < mat.num_rows()-1) UG_LOG(", ");
	}
	UG_LOG("]);\n");
	UG_LOG("(" << mat.num_rows() << " x " << mat.num_cols() << ")\n");

}

template<typename T>
void MemSwap(T &a, T &b)
{
	char c[sizeof(T)];
	memcpy(c, &a, sizeof(T));
	memcpy(&a, &b, sizeof(T));
	memcpy(&b, c, sizeof(T));
}

} // namespace ug

#endif