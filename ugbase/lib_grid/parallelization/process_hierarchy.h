/*
 * Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
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

// NOTE: Classes in this file were originally defined in 'load_balancer.h'
#ifndef __H__UG_process_hierarchy
#define __H__UG_process_hierarchy

#include <vector>
#include "common/util/variant.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl_process_communicator.h"
#endif

namespace ug{

class ProcessHierarchy;
using SPProcessHierarchy = SmartPtr<ProcessHierarchy>;
using ConstSPProcessHierarchy = ConstSmartPtr<ProcessHierarchy>;

/// \addtogroup lib_grid_parallelization_distribution
///	\{

///	Defines how the different levels of a grid shall be distributed across the available processes
/**	Used by LoadBalancer and by different partitioners.*/
class ProcessHierarchy{
	public:
		static SPProcessHierarchy create() {return SPProcessHierarchy(new ProcessHierarchy);}
		~ProcessHierarchy() = default;

	//todo:	add a proc list for more sophisticated hierarchy levels
		void add_hierarchy_level(size_t gridLvl, size_t numProcsPerProc);

		[[nodiscard]] bool empty() const;
		[[nodiscard]] size_t num_hierarchy_levels() const;
		[[nodiscard]] size_t num_global_procs_involved(size_t hierarchyLevel) const;
		[[nodiscard]] size_t grid_base_level(size_t hierarchyLevel) const;
		[[nodiscard]] size_t hierarchy_level_from_grid_level(size_t gridLevel) const;

	/**	Contains all processes which participate on the given hierarchy level,
	 * but only if the local process participates itself. If it doesn't, the
	 * returned communicator is empty.*/
		[[nodiscard]] pcl::ProcessCommunicator global_proc_com(size_t hierarchyLevel) const;

	/**	Contains only processes which are contained in the cluster of the given
	 * hierarchyLevel in which the local process is included.*/
		//pcl::ProcessCommunicator cluster_proc_com(size_t hierarchyLevel);
		[[nodiscard]] const std::vector<int>& cluster_procs(size_t hierarchyLevel) const;

	///	Returns a string which describes the hierarchy layout.
		[[nodiscard]] std::string to_string() const;

	///	allows specification of partitioning hints per hierarchy level.
	/**	\note	A hint is valid for all subsequent hierarchy levels until it is 
	 *			specified again with a different value.
	 *
	 * \note 	Partitioners are free to ignore any partitioning hints.
	 *			For a list of supported hints please check the documentation
	 *			of the respective partitioner.*/
		void add_partition_hint(int hlvl, const std::string& name, const Variant& value);

	///	returns true if the queried partition hint exists and writes it value to 'valOut'
	/**	The method searches hierarchy levels starting from the specified one
	 * down to level 0 and returns the first matching value.*/
		bool partition_hint(Variant& valOut, int hlvl, const std::string& name) const;

	protected:
		using PartitionHintMap = std::map<std::string, Variant>;

		struct HLevelInfo{
			pcl::ProcessCommunicator globalCom;
			std::vector<int> clusterProcs;
			size_t gridLvl;
			size_t numGlobalProcsInUse;
			PartitionHintMap m_partitionHints;
		};

		HLevelInfo& get_hlevel_info(size_t lvl) {return m_levels.at(lvl);}
		[[nodiscard]] const HLevelInfo& get_hlevel_info(size_t lvl) const {return m_levels.at(lvl);}

		void init_cluster_procs(std::vector<int>& clusterProcs,
								size_t hlvl,
								size_t numProcsPerProc);

	private:
		std::vector<HLevelInfo> m_levels;
};

///	\}
	
}//	end of namespace

#endif