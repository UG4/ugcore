// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 22, 2013 (d,m,y)

#ifndef __H__UG__load_balancer__
#define __H__UG__load_balancer__

#include <vector>
#include "lib_grid/multi_grid.h"
#include "lib_grid/tools/partition_map.h"
#include "lib_grid/algorithms/serialization.h"
#include "common/util/table.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl_process_communicator.h"
#endif


namespace ug{

class ProcessHierarchy;
typedef SmartPtr<ProcessHierarchy> SPProcessHierarchy;
typedef ConstSmartPtr<ProcessHierarchy> ConstSPProcessHierarchy;

/// \addtogroup lib_grid_parallelization_distribution
///	\{

///	Defines how the different levels of a grid shall be distributed across the available processes
/**	Used by LoadBalancer and by different partitioners.*/
class ProcessHierarchy{
	public:
		static SPProcessHierarchy create()	{return SPProcessHierarchy(new ProcessHierarchy);}
		~ProcessHierarchy();

	//todo:	add a proc list for more sophisticated hierarchy levels
		void add_hierarchy_level(size_t gridLvl, size_t numProcsPerProc);

		bool empty() const;
		size_t num_hierarchy_levels() const;
		size_t num_global_procs_involved(size_t hierarchyLevel) const;
		size_t grid_base_level(size_t hierarchyLevel) const;
		size_t hierarchy_level_from_grid_level(size_t gridLevel) const;

	/**	Contains all processes which participate on the given hierarchy level,
	 * but only if the local process participates itself. If it doesn't, the
	 * returned communicator is empty.*/
		pcl::ProcessCommunicator global_proc_com(size_t hierarchyLevel) const;

	/**	Contains only processes which are contained in the cluster of the given
	 * hierarchyLevel in which the local process is included.*/
		//pcl::ProcessCommunicator cluster_proc_com(size_t hierarchyLevel);
		const std::vector<int>& cluster_procs(size_t hierarchyLevel) const;

	///	Returns a string which describes the hierarchy layout.
		std::string to_string() const;

	protected:
		struct HLevelInfo{
			pcl::ProcessCommunicator globalCom;
			//pcl::ProcessCommunicator clusterCom;
			std::vector<int> clusterProcs;
			size_t gridLvl;
			size_t numGlobalProcsInUse;
		};

		const HLevelInfo& get_hlevel_info(size_t lvl) const	{return m_levels.at(lvl);}

//		virtual pcl::ProcessCommunicator
//			create_cluster_communicator(size_t hlvl, size_t gridLvl, size_t numProcsPerProc);
		void init_cluster_procs(std::vector<int>& clusterProcs, size_t hlvl,
								size_t numProcsPerProc);

	private:
		std::vector<HLevelInfo>	m_levels;
};


//template <int dim>
//class ConnectionWeights{
//	public:
//		virtual ~ConnectionWeights()	{}
//		virtual void refresh_weights(int baseLevel) = 0;
//		virtual number get_weight(Vertex*) = 0;
//		virtual number get_weight(Edge*) = 0;
//		virtual number get_weight(Face*) = 0;
//};


class IBalanceWeights{
	public:

		virtual ~IBalanceWeights()	{}
		virtual void refresh_weights(int baseLevel)	{};

		virtual number get_weight(Vertex*) = 0;
		virtual number get_weight(Edge*) = 0;
		virtual number get_weight(Face*) = 0;
		virtual number get_weight(Volume*) = 0;
};

typedef SmartPtr<IBalanceWeights>		SPBalanceWeights;


///	Partitioners can be used inside a LoadBalancer or separately to create partition maps
/**	This is the abstract base class for partitioners.*/
class IPartitioner{
	public:
		IPartitioner() :
			m_verbose(true),
			m_clusteredSiblings(true)	{}

		virtual ~IPartitioner()	{}

		virtual void set_next_process_hierarchy(SPProcessHierarchy procHierarchy) = 0;
		virtual void set_balance_weights(SPBalanceWeights balanceWeights) = 0;
//		virtual void set_connection_weights(SmartPtr<ConnectionWeights<dim> > conWeights) = 0;

		virtual ConstSPProcessHierarchy current_process_hierarchy() const = 0;
		virtual ConstSPProcessHierarchy next_process_hierarchy() const = 0;

		virtual bool supports_balance_weights() const = 0;
//		virtual bool supports_connection_weights() const = 0;
		virtual bool supports_repartitioning() const = 0;

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
		virtual void enable_clustered_siblings(bool bEnable)	{m_clusteredSiblings = bEnable;}
		virtual bool clustered_siblings_enabled()				{return m_clusteredSiblings;}
	/**	\} */

	/** The returned distribution quality represents the global quality and thus
	 * is the same for all processes.
	 * You may optionally specify the pointer to a std::vector which will be filled
	 * with the distribution-qualities for each level. By default the pointer is
	 * set to NULL and no level-qualities are thus returned.
	 * If a process doesn't participate on a given level and process, it will write -1
	 * to the corresponding entry in pLvlQualitiesOut.*/
		virtual number estimate_distribution_quality(std::vector<number>* pLvlQualitiesOut = NULL) = 0;

	/**	If the partitioner returns false, no partition-map has been created and
	 * no redistribution should be performed.
	 * Note, that returning false is not a sign of an error - it may even be a
	 * feature of a particular partitioner.*/
		virtual bool partition(size_t baseLvl, size_t elementThreshold) = 0;

		virtual SubsetHandler& get_partitions() = 0;


	///	returns the processes map. Updated during partitioning. may be NULL.
	/**	If NULL is returned, this means that each subset index correspons to a
	 * global proc-rank.*/
		virtual const std::vector<int>* get_process_map() const = 0;

		void set_verbose(bool verbose)	{m_verbose = verbose;}
		bool verbose() const			{return m_verbose;}

	private:
		bool m_verbose;
		bool m_clusteredSiblings;
};

typedef SmartPtr<IPartitioner>		SPPartitioner;


///	A load-balancer redistributes grids using the specified partitioner and process-hierarchy
class LoadBalancer{
	public:
		LoadBalancer();

		virtual ~LoadBalancer();

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
	 * edges which connect elements on different processes.*/
		virtual void rebalance();


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

		void create_quality_record(const char* label);
		void print_quality_records() const;

	private:
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
