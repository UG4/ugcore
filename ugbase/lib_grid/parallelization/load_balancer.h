// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 22, 2013 (d,m,y)

#ifndef __H__UG__load_balancer__
#define __H__UG__load_balancer__

#include <vector>
#include "lib_grid/multi_grid.h"
#include "lib_grid/tools/partition_map.h"
#include "lib_grid/algorithms/serialization.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl_process_communicator.h"
#endif


namespace ug{

class ProcessHierarchy{
	public:
		virtual ~ProcessHierarchy();

	//todo:	add a proc list for more sophisticated hierarchy levels
		void add_hierarchy_level(size_t gridLvl, size_t numProcsPerProc);

		size_t num_hierarchy_levels();
		size_t num_global_procs_involved(size_t hierarchyLevel);
		size_t grid_base_level(size_t hierarchyLevel);
		size_t hierarchy_level_from_grid_level(size_t gridLevel);

	/**	Contains all processes which participate on the given hierarchy level,
	 * but only if the local process participates itself. If it doesn't, the
	 * returned communicator is empty.*/
	//	ProcessCommunicator global_proc_com(size_t hierarchyLevel);

	/**	Contains only processes which are contained in the cluster of the given
	 * hierarchyLevel in which the local process is included.*/
		//pcl::ProcessCommunicator cluster_proc_com(size_t hierarchyLevel);

	protected:
		struct HLevelInfo{
			//ProcessCommunicator globalCom;
			//pcl::ProcessCommunicator clusterCom;
			size_t gridLvl;
			size_t numGlobalProcsInUse;
		};

		const HLevelInfo& get_hlevel_info(size_t lvl) const	{return m_levels.at(lvl);}

//		virtual pcl::ProcessCommunicator
//			create_cluster_communicator(size_t hlvl, size_t gridLvl, size_t numProcsPerProc);

	private:
		std::vector<HLevelInfo>	m_levels;
};

typedef SmartPtr<ProcessHierarchy> SPProcessHierarchy;


template <int dim>
class ConnectionWeights{
	public:
		typedef typename GeomObjBaseTypeByDim<dim>::base_obj_type	elem_type;

		virtual ~ConnectionWeights()	{}
		virtual void set_grid(MultiGrid* mg, Attachment<MathVector<dim> > aPos) = 0;
		virtual void refresh_weights(int baseLevel) = 0;
		virtual number get_weight(elem_type* e1, elem_type* e2) = 0;
	//todo: a get_weight method for parallel environments is required. e.g.
	//		through get_weight(elem_type* e, int side)
};


template <int dim>
class BalanceWeights{
	public:
		typedef typename GeomObjBaseTypeByDim<dim>::base_obj_type	elem_type;

		virtual ~BalanceWeights()	{}
		virtual void set_grid(MultiGrid* mg, Attachment<MathVector<dim> > aPos) = 0;
		virtual void refresh_weights(int baseLevel) = 0;
		virtual number get_weight(elem_type* e) = 0;
};


template <int dim>
class IPartitioner{
	public:
		typedef typename GeomObjBaseTypeByDim<dim>::base_obj_type	elem_t;

		virtual ~IPartitioner()	{}

		virtual void set_grid(MultiGrid* mg, Attachment<MathVector<dim> > aPos) = 0;
		virtual void set_process_hierarchy(SPProcessHierarchy procHierarchy) = 0;
		virtual void set_balance_weights(SmartPtr<BalanceWeights<dim> > balanceWeights) = 0;
		virtual void set_connection_weights(SmartPtr<ConnectionWeights<dim> > conWeights) = 0;

		virtual bool supports_balance_weights() = 0;
		virtual bool supports_connection_weights() = 0;

		virtual void partition(size_t baseLvl) = 0;

		virtual SubsetHandler& get_partitions() = 0;
};


template <int dim>
class LoadBalancer{
	public:
		typedef typename GeomObjBaseTypeByDim<dim>::base_obj_type	elem_t;

		LoadBalancer();

		virtual ~LoadBalancer();

		virtual void set_grid(MultiGrid* mg, Attachment<MathVector<dim> > aPos);

		virtual void enable_vertical_interface_creation(bool enable);

	///	Sets the partitioner which is used to partition the grid into balanced parts.
		virtual void set_partitioner(SmartPtr<IPartitioner<dim> > partitioner);

	///	Sets a callback class which provides the balance weight for a given element
	/**	Balance weights are used to calculate the current balance and to specify
	 * the weight of an element during redistribution.
	 * \note balance weights are only used if the given partitioner supports them.*/
	 	virtual void set_balance_weights(SmartPtr<BalanceWeights<dim> > balanceWeights);

	///	Sets a callback class which provides connection weights between two neighbor elements
	/**	The higher the weight, the higher the likelyhood that the two elements
	 * will reside on the same process after redistribution.
	 * \note connection weights are only used if the given partitioner supports them.*/
	 	virtual void set_connection_weights(SmartPtr<ConnectionWeights<dim> > conWeights);

	///	Inserts a new distribution level on which the grid may be redistributed
	/** Use this method to map a region of levels to a subset of the active processes.
	 * Very useful for hierarchical distribution. Elements between the given level
	 * and the level specified in the next cut are distributed between the given
	 * number of processes only.
	 * The new level has to lie above the level added before.*/
		virtual void add_distribution_level(size_t lvl, size_t numProcsPerProc);

	///	performs load balancing if the balance is too bad or if a distribution level has been reached.
	/**	The balance is calculated using the provieded BalanceWeights class. If
	 * an imbalance is detected (minBalanceWeight / maxBalanceWeight < balanceThreshold)
	 * on a given set of processes, then redistribution will be performed.
	 *
	 * If a distribution level has been reached, e.g. due to refienement, the grid
	 * will be redistributed on that level even if the balance is perfectly fine.
	 *
	 * During redistribution the LoadBalancer tries to distribute elements such that
	 * the sum of balance weights of elements on each process is the same on a
	 * given level. Furthermore it tries to minimize the connection-weights of
	 * edges which connect elements on different processes.*/
		virtual void rebalance(number balanceThreshold);


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

	private:
		typedef Attachment<MathVector<dim> >		APos;
		typedef SmartPtr<IPartitioner<dim> >		SPPartitioner;
		typedef SmartPtr<BalanceWeights<dim> >		SPBalanceWeights;
		typedef SmartPtr<ConnectionWeights<dim> >	SPConnectionWeights;

		MultiGrid*			m_mg;
		APos				m_aPos;
		SPProcessHierarchy	m_processHierarchy;
		SPPartitioner		m_partitioner;
		SPBalanceWeights	m_balanceWeights;
		SPConnectionWeights	m_connectionWeights;
		GridDataSerializationHandler	m_serializer;
		bool m_createVerticalInterfaces;
};

}	// end of namespace

#endif
