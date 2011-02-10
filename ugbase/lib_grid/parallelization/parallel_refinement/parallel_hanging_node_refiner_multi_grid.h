// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 09.02.2011 (m,d,y)

#ifndef __H__UG__parallel_hanging_node_refiner_multi_grid__
#define __H__UG__parallel_hanging_node_refiner_multi_grid__

#include "lib_grid/lg_base.h"
#include "lib_grid/multi_grid.h"
#include "lib_grid/algorithms/refinement/hanging_node_refiner_multi_grid.h"
#include "../distributed_grid.h"
#include "pcl/pcl_communicator.h"

namespace ug
{

/**	The parallel hanging node refiner allows parallel refinement
 * with hanging nodes. Make sure that you initialize it with
 * a valid DistributedGridManager.
 */
class ParallelHangingNodeRefiner_MultiGrid :
	public HangingNodeRefiner_MultiGrid
{
	typedef class HangingNodeRefiner_MultiGrid BaseClass;
	using BaseClass::mark_for_refinement;

	public:
		ParallelHangingNodeRefiner_MultiGrid(IRefinementCallback* refCallback = NULL);

		ParallelHangingNodeRefiner_MultiGrid(
				DistributedGridManager& distGridMgr,
				IRefinementCallback* refCallback = NULL);

		virtual ~ParallelHangingNodeRefiner_MultiGrid();

		void set_distributed_grid_manager(DistributedGridManager& distGridMgr);

	///	all marks are cleared and flags are resetted.
		virtual void clear_marks();

	///	Marks an edge for refinement.
	/**	If interface elements are selected a flag will be checked.*/
		virtual void mark_for_refinement(EdgeBase* e);

	///	Marks a face for refinement.
	/**	If interface elements are selected a flag will be checked.*/
		virtual void mark_for_refinement(Face* f);

	///	Marks a volume for refinement.
	/**	If interface elements are selected a flag will be checked.*/
		virtual void mark_for_refinement(Volume* v);

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
	///	prepares selection and calls the base implementation
	/**	Makes sure that no elements with children are selected,
	 *  Additionally vertices are marked, which have no children but
	 *  associated elements which are marked for refinement.
	 *
	 *  Parallel communication is performed here to notify
	 *  neighboring processes about refinement marks.*/
		virtual void collect_objects_for_refine();

	///	creates required vertices in higher levels.
	/**	Notifies the associated distGridMgr that new elements
	 * may now be created.*/
		virtual void pre_refine();

	/**	Notifies the associated distGridMgr that new elements
	 * have been created.*/
		virtual void post_refine();

	private:
		DistributedGridManager* m_pDistGridMgr;
		pcl::ProcessCommunicator m_procCom;
		pcl::ParallelCommunicator<EdgeLayout> m_intfComEDGE;
		pcl::ParallelCommunicator<FaceLayout> m_intfComFACE;

		bool m_bNewInterfaceEdgesMarked;
		bool m_bNewInterfaceFacesMarked;
		bool m_bNewInterfaceVolumesMarked;
};

}//	end of namespace

#endif
