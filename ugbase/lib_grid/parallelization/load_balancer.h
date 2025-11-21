/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__load_balancer__
#define __H__UG__load_balancer__

#include <vector>
#include "lib_grid/multi_grid.h"
#include "lib_grid/algorithms/serialization.h"
#include "common/util/table.h"
#include "partitioner.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl_process_communicator.h"
#endif


namespace ug{


///	A load-balancer redistributes grids using the specified partitioner and process-hierarchy
class LoadBalancer{
	public:
		LoadBalancer();

		virtual ~LoadBalancer() = default;

		virtual void set_grid(MultiGrid* mg);

		virtual void enable_vertical_interface_creation(bool enable);

	///	Sets the partitioner which is used to partition the grid into balanced parts.
		virtual void set_partitioner(SPPartitioner partitioner);

	///	Sets a callback class which provides the balance weight for a given element
	/**	Balance weights are used to calculate the current balance and to specify
	 * the weight of an element during redistribution.
	 * \note balance weights are only used if the given partitioner supports them.*/
	 	virtual void set_balance_weights(SPBalanceWeights balanceWeights);

	///	Sets a callback class which provides connection weights between two neighbor elements
	/**	The higher the weight, the higher the likelyhood that the two elements
	 * will reside on the same process after redistribution.
	 * \note connection weights are only used if the given partitioner supports them.*/
//	 	virtual void set_connection_weights(SmartPtr<IConnectionWeights> conWeights);

//	///	Inserts a new distribution level on which the grid may be redistributed
//	/** Use this method to map a region of levels to a subset of the active processes.
//	 * Very useful for hierarchical distribution. Elements between the given level
//	 * and the level specified in the next cut are distributed between the given
//	 * number of processes only.
//	 * The new level has to lie above the level added before.*/
//		virtual void add_distribution_level(size_t lvl, size_t numProcsPerProc);

	///	Defines the process hierarchy which will be used during the following calls of rebalance
	 	virtual void set_next_process_hierarchy(SPProcessHierarchy procHierarchy);

	///	If the balance falls below the given threshold, then rebalance will perform redistribution
	/**	Set to 0.9 by default.*/
		virtual void set_balance_threshold(number threshold);

	/**	If distribution on a given level would lead to less elements per process
	 * than the given threshold (in average), then no redistribution will be
	 * performed on that level. Default is 1.*/
		virtual void set_element_threshold(size_t threshold);

//	///	returns the quality of the current distribution
//		virtual number distribution_quality();

	///	performs load balancing if the balance is too bad or if a distribution level has been reached.
	/**	The balance is calculated using the provieded BalanceWeights class. If
	 * an imbalance is detected (minBalanceWeight / maxBalanceWeight < balanceThreshold)
	 * on a given set of processes, then redistribution will be performed.
	 *
	 * If a distribution level has been reached, e.g. due to refinement, the grid
	 * will be redistributed on that level even if the balance is perfectly fine.
	 *
	 * Note that no redistribution is performed on a given level, if the number of
	 * elements on the target process would in average be below the element threshold.
	 *
	 * During redistribution the LoadBalancer tries to distribute elements such that
	 * the sum of balance weights of elements on each process is the same on a
	 * given level. Furthermore it tries to minimize the connection-weights of
	 * edges which connect elements on different processes.
	 *
	 * The method returns false if e.g. problems during partitioning occurred.*/
		virtual bool rebalance();


	/** The returned distribution quality represents the global quality of the elements
	 * of highest dimension and is the same on all processes.
	 * You may optionally specify a pointer to a std::vector which will be filled
	 * with the distribution-qualities for each level. By default the pointer is
	 * set to nullptr and no level-qualities are thus returned.
	 * If a process doesn't participate on a given level, it will write -1
	 * to the corresponding entry in pLvlQualitiesOut.
	 * \{ */
		virtual number estimate_distribution_quality(std::vector<number>* pLvlQualitiesOut);
		number estimate_distribution_quality();
	/** \} */

	///	add serialization callbacks.
	/** Used when the grid is being distributed to pack data associated with grid
	 * objects or associated classes like subset-handlers.
	 * \{ */
		virtual void add_serializer(SPVertexDataSerializer cb)	{m_serializer.add(cb);}
		virtual void add_serializer(SPEdgeDataSerializer cb)	{m_serializer.add(cb);}
		virtual void add_serializer(SPFaceDataSerializer cb)	{m_serializer.add(cb);}
		virtual void add_serializer(SPVolumeDataSerializer cb)	{m_serializer.add(cb);}
		virtual void add_serializer(SPGridDataSerializer cb)	{m_serializer.add(cb);}
	/**	\} */

	///	indicates whether problems occurred during the last rebalancing.
	/**	This can e.g. happen if the partitioner has problems to create a
	 * partition-map which fullfills the given specifications.
	 * \note	if rebalance returns true, still performed rebalancing, even
	 *			if problems occurred. However, the new partition may not
	 *			fullfill all given specifications.*/
		bool problems_occurred();

		void create_quality_record(const char* label);
		void print_quality_records() const;
		void print_last_quality_record() const;

	private:
		template <typename TElem>
		number estimate_distribution_quality_impl(std::vector<number>* pLvlQualitiesOut);

		MultiGrid*			m_mg;
		number				m_balanceThreshold;
		size_t				m_elementThreshold;
		SPProcessHierarchy	m_processHierarchy;
		SPPartitioner		m_partitioner;
		SPBalanceWeights	m_balanceWeights;
//		SPConnectionWeights	m_connectionWeights;
		GridDataSerializationHandler	m_serializer;
		StringStreamTable	m_qualityRecords;
		bool m_createVerticalInterfaces;
};

///	\}
}	// end of namespace

#endif
