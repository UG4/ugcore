øunused
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

#ifndef IVECTOR_VEC_FUNCTIONS_H_
#define IVECTOR_VEC_FUNCTIONS_H_


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// T, IVector functions
template<typename T>
inline void VecAdd(double a1, T &v1, double a2, const IVector &v2)
{
	VecAdd(a1, v1, a2, v2.downcast<T>());
}

template<typename T>
inline void VecAdd(double a1, T &v1, double a2, const IVector &v2, double a3, const IVector &v3)
{
	VecAdd(a1, v1, a2, v2.downcast<T>(), a3, v3.downcast<T>());
}

template<typename T>
inline void VecAdd(double a1, T &v1, double a2, const IVector &v2, double a3, const IVector &v3, double a4, const IVector &v4)
{
	VecAdd(a1, v1, a2, v2.downcast<T>(), a3, v3.downcast<T>(), a4, v4.downcast<T>());
}

template<typename T>
inline double VecProd(T &v1, IVector &v2)
{
	return VecProd(v1, v2.downcast<T>());
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// IVector, IVector functions
void VecAdd(double a1, IVector &v1, double a2, const IVector &v2)
{
	v1.vec_add(a1, a2, v2);
}

void VecAdd(double a1, IVector &v1, double a2, const IVector &v2, double a3, const IVector &v3)
{
	v1.vec_add(a1, a2, v2, a3, v3);
}

void VecAdd(double a1, IVector &v1, double a2, const IVector &v2, double a3, const IVector &v3, double a4, const IVector &v4)
{
	v1.vec_add(a1, a2, v2, a3, v3, a4, v4);
}

double VecProd(const IVector &v1, const IVector &v2)
{
	return v1.vec_prod(v2);
}

double VecNorm2(const IVector &v1)
{
	return v1.norm2();
}

double VecNorm(const IVector &v1)
{
	return v1.norm();
}


#endif