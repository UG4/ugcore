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

#ifndef NUMBER_UTIL_H_
#define NUMBER_UTIL_H_

#include "common/math/ugmath.h"
namespace ug{
extern bool g_bNoNANCheck;
inline bool IsFiniteAndNotTooBig(double d)
{
	const double tooBigValue = 1e30;
	if(d > tooBigValue || d < -tooBigValue || std::isfinite(d) == false)
//	if(std::isfinite(d) == false)
		return false;
	return true;
}

inline bool CloseToZero(double d)
{
	if(d > -1e-100 && d < +1e-100)
		return true;
	else
		return false;
}

//inline bool IsFiniteAndNotTooBig(float d)
//{
//	return IsFiniteAndNotTooBig((double) d);
//}


template <std::size_t N, std::size_t M, typename T>
inline bool IsFiniteAndNotTooBig(const MathMatrix<N, M, T> &m)
{
	for(size_t r=0; r<m.num_rows(); r++)
		for(size_t c=0; c<m.num_cols(); c++)
			if(IsFiniteAndNotTooBig(m(r,c)) == false) return false;
	return true;
}

template <std::size_t N, typename T>
inline bool IsFiniteAndNotTooBig(const MathVector<N, T> &m)
{
	for(size_t r=0; r<m.size(); r++)
			if(IsFiniteAndNotTooBig(m[r]) == false) return false;
	return true;
}



template <size_t TRank, size_t N, typename T>
inline bool IsFiniteAndNotTooBig(const MathTensor<TRank, N, T> &t)
{
	for(size_t i=0; i<t.size(); i++)
		if(IsFiniteAndNotTooBig(t[i]) == false) return false;

	return true;
}


template<typename TData, size_t N>
inline bool IsFiniteAndNotTooBig(const MathTensorX<TData, N> &t)
{
	for(size_t i=0; i<t.size(); i++)
		if(IsFiniteAndNotTooBig(t[i]) == false) return false;

	return true;
}

template<typename TData>
inline bool IsFiniteAndNotTooBig(const std::vector<TData> &t)
{
	for(size_t i=0; i<t.size(); i++)
		if(IsFiniteAndNotTooBig(t[i]) == false) return false;

	return true;
}

template<typename T>
inline bool IsVectorFiniteAndNotTooBig(const T &t)
{
	for(size_t i=0; i<t.size(); i++)
		if(IsFiniteAndNotTooBig(t[i]) == false) return false;

	return true;
}


template<typename T>
inline bool WarnIfIsVectorNanOrTooBig(const T &t, const char *callerName)
{
	if(g_bNoNANCheck) return true;
	size_t iFirst, iCount=0;
	for(size_t i=0; i<t.size(); i++)
		if(IsFiniteAndNotTooBig(t[i]) == false)
		{
			if(iCount == 0) iFirst = i;
			iCount++;
		}

	if(iCount > 0)
	{
		UG_LOG("WARNING in " << callerName << ": " << iFirst << " nan or too big (" << t[iFirst] << "). ");
		if(iCount > 1)
		{ UG_LOG(iCount-1 << " more.");}
		UG_LOG("\n");
		return false;
	}
	return true;
}


template<typename T>
inline bool ThrowIfIsVectorNanOrTooBig(const T &t, const char *callerName)
{
	if(g_bNoNANCheck) return true;
	size_t iFirst, iCount=0;
	for(size_t i=0; i<t.size(); i++)
		if(IsFiniteAndNotTooBig(t[i]) == false)
		{
			if(iCount == 0) iFirst = i;
			iCount++;
		}

	if(iCount > 0)
	{
		UG_THROW("In " << callerName << ": " << iFirst << " nan or too big (" << t[iFirst] << "). "
				<< iCount-1 << " more.");
		return false;
	}

	return true;
}

}


#endif
