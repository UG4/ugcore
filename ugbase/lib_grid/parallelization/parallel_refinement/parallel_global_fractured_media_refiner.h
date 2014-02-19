//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m04 d14

#ifndef __H__LIB_GRID__PARALLEL_GLOBAL_FRACTURED_MEDIA_REFINER__
#define __H__LIB_GRID__PARALLEL_GLOBAL_FRACTURED_MEDIA_REFINER__

#include "../distributed_grid.h"
#include "lib_grid/algorithms/refinement/global_fractured_media_refiner.h"
#include "pcl/pcl_interface_communicator.h"

namespace ug
{

/// \addtogroup lib_grid_parallelization_refinement
/// @{

///	Adds parallel support to GlobalFracturedMediaRefiner
class ParallelGlobalFracturedMediaRefiner : public GlobalFracturedMediaRefiner
{
	public:
		ParallelGlobalFracturedMediaRefiner(DistributedGridManager& distGridMgr);
		virtual ~ParallelGlobalFracturedMediaRefiner();

	protected:
		virtual bool refinement_is_allowed(Vertex* elem);
		virtual bool refinement_is_allowed(EdgeBase* elem);
		virtual bool refinement_is_allowed(Face* elem);
		virtual bool refinement_is_allowed(Volume* elem);
		
		virtual void refinement_step_begins();
		virtual void refinement_step_ends();

		virtual void communicate_marks(BoolMarker& marker);

	protected:
		DistributedGridManager& m_distGridMgr;
		pcl::InterfaceCommunicator<EdgeLayout> m_intfComEDGE;
		pcl::InterfaceCommunicator<FaceLayout> m_intfComFACE;
};

/// @}
}//	end of namespace

////////////////////////////////
//	include implementation
#include "parallel_global_refiner_t_impl.hpp"

#endif
