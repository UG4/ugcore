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

#ifndef __H__UG_partitioner
#define __H__UG_partitioner

#include "process_hierarchy.h"
#include "lib_grid/multi_grid.h"
#include "lib_grid/tools/subset_handler_grid.h"

namespace ug {

/// \addtogroup lib_grid_parallelization_distribution
///	\{

class IPartitioner;

class IBalanceWeights{
	public:

		virtual ~IBalanceWeights() = default;
		virtual void refresh_weights(int baseLevel)	{};

		virtual number get_weight(Vertex*) {return 1.;}
		virtual number get_weight(Edge*) {return 1.;}
		virtual number get_weight(Face*) {return 1.;}
		virtual number get_weight(Volume*) {return 1.;}

		virtual number get_refined_weight(Vertex* e) {return get_weight(e);}
		virtual number get_refined_weight(Edge* e) {return 2. * get_weight(e);}
		virtual number get_refined_weight(Face* e) {return 4. * get_weight(e);}///< todo: use a more sophisticated implementation
		virtual number get_refined_weight(Volume* e) {return 8. * get_weight(e);}///< todo: use a more sophisticated implementation
		
		virtual bool has_level_offsets()		{return false;}

	///	Relative indicator in which level the specified elements should be partitioned.
	/** If this method returns true, one should use get_refined_weight instead of get_weight.
	 * \{ */
		virtual bool consider_in_level_above(Vertex*) {return false;}
		virtual bool consider_in_level_above(Edge*) {return false;}
		virtual bool consider_in_level_above(Face*) {return false;}
		virtual bool consider_in_level_above(Volume*) {return false;}
	/** \} */
};

using SPBalanceWeights = SmartPtr<IBalanceWeights>;


/**
 * Interface for the definition of weights for connections between elements.
 * These weights introduce a customizable penalty for any side if the border
 * of a partition includes the side.
 * This can be used to minimize communication cost e.g. when using the ParMetis
 * partitioner.
 * It can also be useful if you want to make sure that certain subsets are not
 * contained in any partition boundary or if you want to enforce that a collection
 * of elements be in one single partition.
 */
class ICommunicationWeights
{
	public:
		virtual ~ICommunicationWeights() = default;

		/**
		 * Get the weight of a specific connection.
		 * @param conn pointer to the connection in question
		 * @return weight for the side
		 */
		virtual number get_weight(GridObject* conn) = 0;

		/**
		 * Whether the given connection is to be assigned another weight.
		 * If true, this weight can be obtained by get_weight(),
		 * if false, no such weight must be requested.
		 * @param conn pointer to the connection in question
		 * @return true iff get_weight() can provide a new weight for the connection
		 */
		virtual bool reweigh(GridObject* conn) {return false;}
};

using SPCommunicationWeights = SmartPtr<ICommunicationWeights>;


///	allows to pre-process data before partitioning starts
/**	If supported by a partitioner, a pre-processor is called before partitioning
 * starts. This may e.g. be useful to provide an alternative coordinate set.
 */
class IPartitionPreProcessor{
	public:
		virtual ~IPartitionPreProcessor() = default;
		
		virtual void partitioning_starts (	MultiGrid* mg,
		                                 	IPartitioner* partitioner) = 0;
		virtual void partitioning_done (	MultiGrid* mg,
		                                	IPartitioner* partitioner) = 0;
};

using SPPartitionPreProcessor = SmartPtr<IPartitionPreProcessor>;



///	allows to post-process partitions
/**	If supported by a partitioner, a post-processor is called for each partitioned level
 * to allow for the adjustment of partitions following some specific rules.
 *
 * 'init_post_processing' is called before 'post_process' is called for the first time.
 * A specialization can e.g. attach required attachments during this routine.
 * 
 * 'post_process' is called each time partitioning is done for a hierarchy-level.
 * It will be called with the multigrid level-index on which the partitioning was performed.
 * Note that this means that post_process will only be called for some MultiGrid levels.
 *
 * When partitioning is completed, 'partitioning_done' will be called. Use this method
 * to detach any attachments that were attached during 'init_post_processing'*/
class IPartitionPostProcessor{
	public:
		virtual ~IPartitionPostProcessor()	= default;
		
		virtual void init_post_processing(	MultiGrid* mg,
		                                  	SubsetHandler* partitions) = 0;
		virtual void post_process(int lvl) = 0;
		virtual void partitioning_done() = 0;
};

using SPPartitionPostProcessor = SmartPtr<IPartitionPostProcessor>;



///	Partitioners can be used inside a LoadBalancer or separately to create partition maps
/**	This is the abstract base class for partitioners.*/
class IPartitioner{
	public:
		IPartitioner() :
			m_problemsOccurred(false),
			m_verbose(true),
			m_clusteredSiblings(true)	{}

		virtual ~IPartitioner()	= default;

		virtual void set_next_process_hierarchy(SPProcessHierarchy procHierarchy) = 0;

		virtual void set_balance_weights(SPBalanceWeights balanceWeights) {
			UG_COND_THROW(!supports_balance_weights(),
			              "This partitioner does not support balance weights.");
		};

		virtual void set_communication_weights(SPCommunicationWeights commWeights)	{
			UG_COND_THROW(!supports_communication_weights(),
			              "This partitioner does not support communication weights.");
		};

		virtual void set_partition_pre_processor(SPPartitionPreProcessor) {
			UG_THROW("Partition-Pre-Processing is currently not supported by the chosen partitioner.");
		}

		virtual void set_partition_post_processor(SPPartitionPostProcessor){
			UG_THROW("Partition-Post-Processing is currently not supported by the chosen partitioner.");
		}

		[[nodiscard]] virtual ConstSPProcessHierarchy current_process_hierarchy() const = 0;
		[[nodiscard]] virtual ConstSPProcessHierarchy next_process_hierarchy() const = 0;

		[[nodiscard]] virtual bool supports_balance_weights() const {return false;}
		[[nodiscard]] virtual bool supports_communication_weights() const {return false;}
		[[nodiscard]] virtual bool supports_repartitioning() const {return false;}

	/**	clustered siblings help to ensure that all vertices which are connected to
	 * a constrained vertex through are on the same process as the constrained vertex.
	 * If only refinement is performed, it would be sufficient to only cluster
	 * constrained siblings. However, coarsening would be rather complicated in that
	 * case, since it is rather complicated to introduce constrained sibling elements if a
	 * previously unconstrained sibling is not located on the same process...
	 *
	 * If you only perform global refinement, you may safely disable clustered siblings.
	 * The distribution quality is most likely better in this case.
	 * \note	enabled by default
	 * \{ */
		virtual void enable_clustered_siblings(bool bEnable) {m_clusteredSiblings = bEnable;}
		virtual bool clustered_siblings_enabled() {return m_clusteredSiblings;}
	/**	\} */


	/**	If the partitioner returns false, no partition-map has been created and
	 * no redistribution should be performed.
	 * Note, that returning false is not a sign of an error - it may even be a
	 * feature of a particular partitioner.*/
		virtual bool partition(size_t baseLvl, size_t elementThreshold) = 0;

		virtual SubsetHandler& get_partitions() = 0;

		virtual void set_subset_handler(SmartPtr<SubsetHandler>) {
			UG_LOG("IPartitioner::set_subset_handler not implemented");
		};
		virtual void enable_static_partitioning(bool) {
			UG_LOG("IPartitioner::enable_static_partitioning not implemented");
		};

	///	returns the processes map. Updated during partitioning. may be nullptr.
	/**	If nullptr is returned, this means that each subset index corresponds to a
	 * global proc-rank.*/
		[[nodiscard]] virtual const std::vector<int>* get_process_map() const = 0;

	///	indicates whether problems occurred during the last partitioning
	/**	\note	if partition(...) returns true, the partition map is valid,
	 *			even if problems occurred. It may however not be optimal.*/
		virtual bool problems_occurred() {return m_problemsOccurred;}

		void set_verbose(bool verbose) {m_verbose = verbose;}
		[[nodiscard]] bool verbose() const {return m_verbose;}

	protected:
		bool m_problemsOccurred;

	private:
		bool m_verbose;
		bool m_clusteredSiblings;
};

using SPPartitioner = SmartPtr<IPartitioner>;

///	\}

}//	end of namespace

#endif