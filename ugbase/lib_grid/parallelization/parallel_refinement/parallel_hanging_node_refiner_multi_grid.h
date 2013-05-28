// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 08.03.2012 (m,d,y)

#ifndef __H__UG__parallel_hanging_node_refiner_multi_grid__
#define __H__UG__parallel_hanging_node_refiner_multi_grid__


#include "lib_grid/algorithms/refinement/hanging_node_refiner_multi_grid.h"
#include "../distributed_grid.h"
#include "pcl/pcl_interface_communicator.h"

namespace ug
{

/// \addtogroup lib_grid_parallelization_refinement
/// @{

/**	This is a template class that allows to use a refiner in a parallel
 * environment.
 * Make sure that you initialize it with a valid DistributedGridManager.
 */
class ParallelHangingNodeRefiner_MultiGrid : public HangingNodeRefiner_MultiGrid
{
	public:
		typedef HangingNodeRefiner_MultiGrid BaseClass;
		using BaseClass::mark;
		using BaseClass::copy_marks_to_vmasters;
		using BaseClass::copy_marks_to_vslaves;

		ParallelHangingNodeRefiner_MultiGrid(IRefinementCallback* refCallback = NULL);

		ParallelHangingNodeRefiner_MultiGrid(
				DistributedGridManager& distGridMgr,
				IRefinementCallback* refCallback = NULL);

		virtual ~ParallelHangingNodeRefiner_MultiGrid();

		void set_distributed_grid_manager(DistributedGridManager& distGridMgr);

	/**	If not all processes are involved in refinement,
	 *	one can set the involved processes here. By default
	 *	all processes are involved.*/
		void set_involved_processes(pcl::ProcessCommunicator com);

	protected:
	///	a callback that allows to deny refinement of special vertices
		virtual bool refinement_is_allowed(VertexBase* elem);
	///	a callback that allows to deny refinement of special edges
		virtual bool refinement_is_allowed(EdgeBase* elem);
	///	a callback that allows to deny refinement of special faces
		virtual bool refinement_is_allowed(Face* elem);
	///	a callback that allows to deny refinement of special volumes
		virtual bool refinement_is_allowed(Volume* elem);

		virtual bool continue_collect_objects_for_refine(bool continueRequired);

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

	/**	Notifies the associated distGridMgr that elements will be erased*/
		virtual void pre_coarsen();

	/**	Notifies the associated distGridMgr that elements have been erased.*/
		virtual void post_coarsen();

	///	copies the current marks in the ref-mark-selector from v-slaves to v-masters
	/**	\{ */
		template <class TElem, class TIntfcCom>
		void copy_marks_to_vmasters(TIntfcCom& com);
	/** \} */

	///	copies the current marks in the ref-mark-selector from v-slaves to v-masters
	/**	\{ */
		template <class TElem, class TIntfcCom>
		void copy_marks_to_vslaves(TIntfcCom& com);
	/** \} */


	///	allows to check whether a distributed grid contains edges
		virtual bool contains_edges();

	///	allows to check whether a distributed grid contains faces
		virtual bool contains_faces();

	///	allows to check whether a distributed grid contains volumes
		virtual bool contains_volumes();

		virtual void broadcast_marks_horizontally(bool vertices, bool edges, bool faces);

		virtual void copy_marks_to_vmasters(bool vertices, bool edges,
											bool faces, bool volumes);

		virtual void copy_marks_to_vslaves(bool vertices, bool edges,
										   bool faces, bool volumes);

		virtual bool one_proc_true(bool localProcTrue);

	private:
		DistributedGridManager* m_pDistGridMgr;
		MultiGrid*				m_pMG;
		pcl::ProcessCommunicator m_procCom;
		pcl::InterfaceCommunicator<VertexLayout> m_intfComVRT;
		pcl::InterfaceCommunicator<EdgeLayout> m_intfComEDGE;
		pcl::InterfaceCommunicator<FaceLayout> m_intfComFACE;
		pcl::InterfaceCommunicator<VolumeLayout> m_intfComVOL;
};

/// @}

}//	end of namespace

////////////////////////////////
//	include implementation

#endif
