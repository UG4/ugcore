/**
 * \file amg_debug_helper.h
 *
 * \author Martin Rupp
 *
 * \date 10.08.2010
 *
 * class declaration for amg
 *
 * Goethe-Center for Scientific Computing 2010.
 */

#ifndef __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_DEBUG_HELPER_H__
#define __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_DEBUG_HELPER_H__

#include "common/common.h" // MathVector<3>
#include "common/math/ugmath.h"

namespace ug{

inline void WritePositionsToStream(std::ostream &out, const std::vector<MathVector<3> > &positions, int dimension)
{
	out << dimension << std::endl;
	out << positions.size() << std::endl;
	for(size_t i=0; i< positions.size() ; i++)
	{
		out << positions[i][0] << " " << positions[i][1];
		if(dimension == 3)
			out << " " << positions[i][2];
		out << std::endl;
	}
}

/**
 * \brief helper struct for debug output of matrices in AMG
 * This class knows the position to each node on each level of
 * the algebraic grid.
 * \sa WriteAMGGraphToFile, AMGWriteToFile, AMGWriteToFilePS
 */
struct cAMG_helper
{
	cAMG_helper() : dimension(0), parentIndex(NULL) {	}

	//! returns to the node i on level 'level' the corresponding node on the level 0
	int GetOriginalIndex(int level, int i) const
	{
		return GetOriginalIndex(level, 0, i);
	}

	//! returns to the node i on level 'level' the corresponding node on the level 0
	int GetOriginalIndex(int level, int onlevel, int i) const
	{
		while(level != onlevel)
		{
			if(i >= (int)(*parentIndex)[level].size())
				return 0;
			i = (*parentIndex)[level--][i];
		}
		return i;
	}

	//! returns to the node i on level 'level' the position
	MathVector<3> GetPosForIndexAtLevel(int level, int i) const
	{
		return positions[level][i];
	}


	bool has_positions() const
	{
		return dimension != 0;
	}

	stdvector<stdvector<MathVector<3> > > positions; ///< positions on the AMG grid 0
	int dimension;					///< dimension (2 or 3)
	stdvector< stdvector<int> > *parentIndex;				///< (*parentIndex)[L][i] is the index of i on level L-1
};


} // namespace ug
#endif // __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_DEBUG_HELPER_H__
