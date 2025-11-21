/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef BLAS_VEC_INTERFACE_H_
#define BLAS_VEC_INTERFACE_H_


// specializations for double
//! calculates dest = alpha1*v1. for doubles
inline void VecScaleAssign(double &dest, double alpha1, const double &v1)
{
	dest = alpha1*v1;
}
// and so on.


//! calculates dest = alpha1*v1
//
template<typename vector_t>
inline void VecScaleAssign(vector_t &dest, double alpha1, const vector_t &v1)
{
	for(size_t i=0; i<dest.size(); i++)
		VecScaleAssign(dest[i], alpha1, v1[i]);
}

//! calculates dest = v1
//
template<typename vector_t>
inline void VecAssign(vector_t &dest, const vector_t &v1)
{
	for(size_t i=0; i<dest.size(); i++)
	    dest[i]=v1[i];
}


//! calculates dest = alpha1*v1 + alpha2*v2
template<typename vector_t>
inline void VecScaleAdd(vector_t &dest, double alpha1, const vector_t &v1, double alpha2, const vector_t &v2);

//! calculates dest = alpha1*v1 + alpha2*v2 + alpha3*v3
template<typename vector_t>
inline void VecScaleAdd(vector_t &dest, double alpha1, const vector_t &v1, double alpha2, const vector_t &v2, double alpha3, const vector_t &v3);


//! calculates dest = dest + v
template<typename vector_t>
inline void VecAdd(vector_t &dest, vector_t &v);

//! calculates dest = dest - v
template<typename vector_t>
inline void VecSubstract(vector_t &dest, vector_t &v);


// VecProd

//! calculates s += scal<a, b>
template<typename vector_t>
inline void VecProdAdd(const vector_t &a, const vector_t &b, double &sum)
{
	for(size_t i=0; i<a.size(); i++) VecProdAdd(a[i], b[i], sum);
}

//! returns scal<a, b>
template<typename vector_t>
inline double VecProd(const vector_t &a, const vector_t &b)
{
	double sum=0;
	VecProdAdd(a, b, sum);
	return sum;
}

// VecNorm

//! calculates s += norm_2^2(a)
template<typename vector_t>
inline void VecNormSquaredAdd(const vector_t &a, const vector_t &b, double &sum)
{
	for(int i=0; i<a.size(); i++) VecNormSquaredAdd(a[i], sum);
}

//! returns norm_2^2(a)
template<typename vector_t>
inline double VecNormSquared(const vector_t &a, const vector_t &b)
{
	double sum=0;
	VecNormSquaredAdd(a, sum);
	return sum;
}

#endif