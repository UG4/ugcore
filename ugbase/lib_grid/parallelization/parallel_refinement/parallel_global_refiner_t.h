//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m04 d14

#ifndef __H__LIB_GRID__PARALLEL_GLOBAL_REFINER_T__
#define __H__LIB_GRID__PARALLEL_GLOBAL_REFINER_T__

#include "../distributed_grid.h"

namespace ug
{

/// \addtogroup lib_grid_parallelization_refinement
/// @{

///	Adds parallel support to a global refiner
/**	\todo	This class would benefit from some comfort methods and a better
 * 			documentation.
 */
template <class TRefiner>
class TParallelGlobalRefiner : public TRefiner
{
	public:
		TParallelGlobalRefiner(DistributedGridManager& distGridMgr);
		virtual ~TParallelGlobalRefiner();

	protected:
		virtual bool refinement_is_allowed(Vertex* elem);
		virtual bool refinement_is_allowed(EdgeBase* elem);
		virtual bool refinement_is_allowed(Face* elem);
		virtual bool refinement_is_allowed(Volume* elem);
		
		virtual void refinement_step_begins();
		virtual void refinement_step_ends();

	protected:
		DistributedGridManager& m_distGridMgr;
};

/// @}
}//	end of namespace

////////////////////////////////
//	include implementation
#include "parallel_global_refiner_t_impl.hpp"

#endif
