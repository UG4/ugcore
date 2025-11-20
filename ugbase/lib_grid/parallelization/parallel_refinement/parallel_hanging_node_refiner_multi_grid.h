/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG__parallel_hanging_node_refiner_multi_grid__
#define __H__UG__parallel_hanging_node_refiner_multi_grid__


#include "lib_grid/refinement/hanging_node_refiner_multi_grid.h"
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
		using BaseClass = HangingNodeRefiner_MultiGrid;
		using BaseClass::mark;
		using BaseClass::copy_marks_to_vmasters;
		using BaseClass::copy_marks_to_vslaves;

		ParallelHangingNodeRefiner_MultiGrid(SPRefinementProjector projector = nullptr);

		ParallelHangingNodeRefiner_MultiGrid(
				DistributedGridManager& distGridMgr,
				SPRefinementProjector projector = nullptr);

		~ParallelHangingNodeRefiner_MultiGrid() override;

		void set_distributed_grid_manager(DistributedGridManager& distGridMgr);

	/**	If not all processes are involved in refinement,
	 *	one can set the involved processes here. By default
	 *	all processes are involved.*/
		void set_involved_processes(pcl::ProcessCommunicator com);

	protected:
	///	a callback that allows to deny refinement of special vertices
		bool refinement_is_allowed(Vertex* elem) override;
	///	a callback that allows to deny refinement of special edges
		bool refinement_is_allowed(Edge* elem) override;
	///	a callback that allows to deny refinement of special faces
		bool refinement_is_allowed(Face* elem) override;
	///	a callback that allows to deny refinement of special volumes
		bool refinement_is_allowed(Volume* elem) override;

		bool continue_collect_objects_for_refine(bool continueRequired) override;

	///	distributes hnode marks
	/**	Calls the base implementation to assign hnode marks and afterwards
	 * distributes them amongst neighbor processes.*/
		void assign_hnode_marks() override;

	///	creates required vertices in higher levels.
	/**	Notifies the associated distGridMgr that new elements
	 * may now be created.*/
		void pre_refine() override;

	/**	Notifies the associated distGridMgr that new elements
	 * have been created.*/
		void post_refine() override;

	/**	Notifies the associated distGridMgr that elements will be erased*/
		void pre_coarsen() override;

	/**	Notifies the associated distGridMgr that elements have been erased.*/
		void post_coarsen() override;

	///	copies the current marks in the ref-mark-selector from v-slaves to v-masters
	/**	\{ */
		template <typename TElem, typename TIntfcCom>
		void copy_marks_to_vmasters(TIntfcCom& com);
	/** \} */

	///	copies the current marks in the ref-mark-selector from v-slaves to v-masters
	/**	\{ */
		template <typename TElem, typename TIntfcCom>
		void copy_marks_to_vslaves(TIntfcCom& com);
	/** \} */


	///	allows to check whether a distributed grid contains edges
		bool contains_edges() override;

	///	allows to check whether a distributed grid contains faces
		bool contains_faces() override;

	///	allows to check whether a distributed grid contains volumes
		bool contains_volumes() override;

		void broadcast_marks_horizontally(bool vertices, bool edges, bool faces,
		                                  bool allowDeselection = false) override;

		void broadcast_marks_vertically(bool vertices, bool edges,
		                                bool faces, bool volumes,
		                                bool allowDeselection = false) override;

		void copy_marks_to_vmasters(bool vertices, bool edges,
		                            bool faces, bool volumes) override;

		void copy_marks_to_vslaves(bool vertices, bool edges,
		                           bool faces, bool volumes) override;

		bool one_proc_true(bool localProcTrue) override;

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
