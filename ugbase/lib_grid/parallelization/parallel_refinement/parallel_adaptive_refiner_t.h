// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 09.02.2011 (m,d,y)

#ifndef __H__UG__parallel_adaptive_refiner_t__
#define __H__UG__parallel_adaptive_refiner_t__

#include "lib_grid/lg_base.h"
#include "lib_grid/multi_grid.h"
#include "lib_grid/algorithms/refinement/hanging_node_refiner_multi_grid.h"
#include "../distributed_grid.h"
#include "pcl/pcl_communicator.h"

namespace ug
{

/// \addtogroup lib_grid_parallelization_refinement
/// @{

/**	This is a template class that allows to use a refiner in a parallel
 * environment.
 * Make sure that you initialize it with a valid DistributedGridManager.
 *
 * \todo	This refiner currently only works for multi-grids.
 */
template <class TRefiner>
class TParallelAdaptiveRefiner:
	public TRefiner
{
	typedef TRefiner BaseClass;
	using BaseClass::mark;

	public:
		TParallelAdaptiveRefiner(IRefinementCallback* refCallback = NULL);

		TParallelAdaptiveRefiner(
				DistributedGridManager& distGridMgr,
				IRefinementCallback* refCallback = NULL);

		virtual ~TParallelAdaptiveRefiner();

		void set_distributed_grid_manager(DistributedGridManager& distGridMgr);

	///	all marks are cleared and flags are resetted.
		virtual void clear_marks();

	///	Marks an vertex for refinement.
		virtual bool mark(VertexBase* v, RefinementMark refMark = RM_REGULAR);

	///	Marks an edge for refinement.
	/**	If interface elements are selected a flag will be checked.*/
		virtual bool mark(EdgeBase* e, RefinementMark refMark = RM_REGULAR);

	///	Marks a face for refinement.
	/**	If interface elements are selected a flag will be checked.*/
		virtual bool mark(Face* f, RefinementMark refMark = RM_REGULAR);

	///	Marks a volume for refinement.
	/**	If interface elements are selected a flag will be checked.*/
		virtual bool mark(Volume* v, RefinementMark refMark = RM_REGULAR);

	/**	If not all processes are involved in refinement,
	 *	one can set the involved processes here. By default
	 *	all processes are involved.*/
		void set_involved_processes(pcl::ProcessCommunicator com);

	///	performs parallel refinement
	/**	Checks that everything was initialized correctly and calls the
	 * base implementation.
	 * Throws an instance of UGError if something went wrong.
	 *
	 * Parallelization is mainly performed in collect_objects_for_refine.*/
		virtual void refine();

	protected:
	///	a callback that allows to deny refinement of special vertices
		virtual bool refinement_is_allowed(VertexBase* elem);
	///	a callback that allows to deny refinement of special edges
		virtual bool refinement_is_allowed(EdgeBase* elem);
	///	a callback that allows to deny refinement of special faces
		virtual bool refinement_is_allowed(Face* elem);
	///	a callback that allows to deny refinement of special volumes
		virtual bool refinement_is_allowed(Volume* elem);

	///	prepares selection and calls the base implementation
	/**	Makes sure that no elements with children are selected,
	 *  Additionally vertices are marked, which have no children but
	 *  associated elements which are marked for refinement.
	 *
	 *  Parallel communication is performed here to notify
	 *  neighboring processes about refinement marks.*/
		virtual void collect_objects_for_refine();

	///	distributes hnode marks
	/**	Calls the base implementation to assign hnode marks and afterwards
	 * distributes them amongst neighbor processes.*/
		virtual void assign_hnode_marks();

	///	creates required vertices in higher levels.
	/**	Notifies the associated distGridMgr that new elements
	 * may now be created.*/
		virtual void pre_refine();

	/**	Notifies the associated distGridMgr that new elements
	 * have been created.*/
		virtual void post_refine();

	private:
		DistributedGridManager* m_pDistGridMgr;
		MultiGrid*				m_pMG;
		pcl::ProcessCommunicator m_procCom;
		pcl::ParallelCommunicator<VertexLayout> m_intfComVRT;
		pcl::ParallelCommunicator<EdgeLayout> m_intfComEDGE;
		pcl::ParallelCommunicator<FaceLayout> m_intfComFACE;

		bool m_bNewInterfaceVerticesMarked;
		bool m_bNewInterfaceEdgesMarked;
		bool m_bNewInterfaceFacesMarked;
		bool m_bNewInterfaceVolumesMarked;
};

/// @}

}//	end of namespace

////////////////////////////////
//	include implementation
#include "parallel_adaptive_refiner_t_impl.hpp"

#endif
