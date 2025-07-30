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

#ifndef VECTOR_UTIL_H_
#define VECTOR_UTIL_H_

#include "lib_algebra/small_algebra/blocks.h"
#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#endif

namespace ug{

template<typename TBlock>
bool BlockVectorFiniteAndNotTooBig(TBlock &v, double tooBigValue=1e24)
{
	for(size_t j=0; j< GetSize(v); j++)
	{
		double d = BlockRef(v, j);
		if(d > tooBigValue || d < -tooBigValue || std::isfinite(d) == false)
				return false;
	}
	return true;
}

template<typename TBlock>
bool BlockMatrixFiniteAndNotTooBig(TBlock &m, double tooBigValue=1e24)
{
	for(size_t r=0; r<GetRows(m); r++)
		for(size_t c=0; c< GetCols(m); c++)
		{
			double d = BlockRef(m, r, c);
			if(d > tooBigValue || d < -tooBigValue || std::isfinite(d) == false)
					return false;
		}

	return true;
}


template<typename TVector>
bool IsFiniteAndNotTooBig(const TVector &v, double tooBigValue=1e24)
{
	for(size_t i=0; i<v.size(); i++)
		if(BlockVectorFiniteAndNotTooBig(v[i], tooBigValue) == false)
			return false;
	return true;
}

/*template<typename TVector>
bool IsFinite(const TVector &v)
{
	for(size_t i=0; i<v.size(); i++)
	{
		for(size_t j=0; j< GetSize(v[i]); j++)
		{
			double d = BlockRef(v[i], j);
			if(isfinite(d) == false)
				return false;
		}
	}
	return true;
}*/

#ifdef UG_PARALLEL
template<typename TVector>
bool IsFiniteAndNotTooBig(const ParallelVector<TVector> &v)
{
	return AllProcsTrue(IsFiniteAndNotTooBig((TVector&)v), v.layouts()->proc_comm());
}
template<typename TVector>
bool IsFinite(const ParallelVector<TVector> &v)
{
	return AllProcsTrue(IsFiniteAndNotTooBig((TVector&)v), v.layouts()->proc_comm());
}
#endif

}
#endif /* VECTOR_UTIL_H_ */
