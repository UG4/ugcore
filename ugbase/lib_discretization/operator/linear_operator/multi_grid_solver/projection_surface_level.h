/*
 * projection_surface_level.h
 *
 *  Created on: 01.10.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__PROJECTION_SURFACE_LEVEL__
#define __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__PROJECTION_SURFACE_LEVEL__

// library intern headers
#include "lib_discretization/function_spaces/grid_function.h"

namespace ug{

template <typename TDoFDistribution, typename TVector>
struct ProjectionSurfaceLevel
{
	static bool surface_to_level(	TVector& levelFunction,
	                             	const IDoFDistribution<TDoFDistribution>& levelDoFDistribution,
									const TVector& surfaceFunction,
									const IDoFDistribution<TDoFDistribution>& surfaceDoFDistribution)
	{
	//  check if DoFDistributions are equal
	//	(i.e. full refinement and surfaceGrid == top of levelGrids)
		if(&levelDoFDistribution != &surfaceDoFDistribution)
		{
			UG_LOG("ERROR in ProjectionSurfaceLevel::surface_to_level(): "
						"Currently only implemented for full refinement\n");
			return false;
		}

	//  copy values
		levelFunction = surfaceFunction;
		return true;
	}

	static bool level_to_surface(	TVector& surfaceFunction,
	                             	const IDoFDistribution<TDoFDistribution>& surfaceDoFDistribution,
									const TVector& levelFunction,
									const IDoFDistribution<TDoFDistribution>& levelDoFDistribution)
	{
	//  check if DoFDistributions are equal
	//	(i.e. full refinement and surfaceGrid == top of levelGrids)
		if(&levelDoFDistribution != &surfaceDoFDistribution)
		{
			UG_LOG("ERROR in ProjectionSurfaceLevel::surface_to_level(): "
						"Currently only implemented for full refinement\n");
			return false;
		}

	//  copy values
		surfaceFunction = levelFunction;
		return true;
	}
};

} // end namespace ug
#endif
