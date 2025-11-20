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

#ifndef VEC_VEC_FUNCTIONS_H_
#define VEC_VEC_FUNCTIONS_H_
#include "common/error.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// T, T functions
template<typename T>
inline void VecAdd(double a1, T &v1, double a2, const T &v2)
{
	ASSERT_EQUAL(v1.size(), v2.size());
	for(size_t i=0; i<v1.size(); i++)
		v1[i] = a1*v1[i] + a2*v2[i];
}

template<typename T>
inline void VecAdd(double a1, T &v1, double a2, const T &v2, double a3, const T &v3)
{
	ASSERT_EQUAL_3(v1.size(), v2.size(), v3.size());
	for(size_t i=0; i<v1.size(); i++)
		v1[i] = a1*v1[i] + a2*v2[i] + a3*v3[i];
}

template<typename T>
inline void VecAdd(double a1, T &v1, double a2, const T &v2, double a3, const T &v3, double a4, const T &v4)
{
	ASSERT_EQUAL_4(v1.size(), v2.size(), v3.size(), v4.size());
	for(size_t i=0; i<v1.size(); i++)
		v1[i] = a1*v1[i] + a2*v2[i] + a3*v3[i] + a4*v4[i];
}

template<typename T>
inline double VecProd(T &v1, T &v2)
{
	ASSERT_EQUAL(v1.size(), v2.size());
	double s=0;
	for(size_t i=0; i<v1.size(); i++)
		s += v1[i]*v2[i];
	return s;
}

template<typename T>
inline double VecNorm2(T &v)
{
	double s=0;
	for(size_t i=0; i<v.size(); i++)
		s += v[i]*v[i];
	return s;
}


#endif