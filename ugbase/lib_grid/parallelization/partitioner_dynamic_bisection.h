// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Nov 4, 2013

#ifndef __H__UG__partitioner_dynamic_biscection__
#define __H__UG__partitioner_dynamic_biscection__

#include <vector>
#include "parallel_grid_layout.h"
#include "load_balancer.h"
#include "pcl/pcl_interface_communicator.h"

namespace ug{

/// \addtogroup lib_grid_parallelization_distribution
///	\{

///	Parallel bisection partitioner
/**	The partitioner can be used inside a LoadBalancer or separately. It can
 * operate on serial and parallel multigrids.
 *
 * The following partition hints provided by a process hierarchy are supported:
 *	- 'enableXCuts' (bool)
 *	- 'enableYCuts' (bool)
 *	- 'enableZCuts' (bool)
 */
template <class TElem, int dim>
class Partitioner_DynamicBisection : public IPartitioner{
	public:
		typedef IPartitioner	 			base_class;
		typedef TElem						elem_t;
		typedef MathVector<dim>				vector_t;
		typedef Attachment<vector_t>		apos_t;
		typedef Grid::VertexAttachmentAccessor<apos_t>	aapos_t;
		typedef typename GridLayoutMap::Types<elem_t>::Layout::LevelLayout	layout_t;

		Partitioner_DynamicBisection();
		virtual ~Partitioner_DynamicBisection();

		void set_grid(MultiGrid* mg, Attachment<MathVector<dim> > aPos);

	///	allows to optionally specify a subset-handler on which the balancer shall operate
		void set_subset_handler(SmartPtr<SubsetHandler> sh);

	///	sets the tolerance threshold. 1: no tolerance, 0: full tolerance.
	/**	the tolerance is defaulted to 0.99*/
		void set_tolerance(number tol)	{m_tolerance = tol;}

	///	enables automatical determination of split-axis by longest geometry extension
	/** disabled by default*/
		void enable_longest_split_axis(bool enable)	{m_longestSplitAxisEnabled = enable;}

	///	enable or disable splits along a certain axis.
	/**	By default, splits along all axis are enabled.*/
		void enable_split_axis(int axis, bool enable);

	///	sets the axis with which bisection is started.
	/** 0 by default.
	 * \note	the start-split-axis is only considered if longest-split-axis is disabled.*/
		void set_start_split_axis(int axis)			{m_startSplitAxis = axis;}

	///	the maximum number of iterations performed to find a good split plane
		int num_split_improvement_iterations() const	{return m_splitImproveIterations;}
		void set_num_split_improvement_iterations(int num)	{m_splitImproveIterations = num;}
		
		virtual void set_next_process_hierarchy(SPProcessHierarchy procHierarchy);
		virtual void set_balance_weights(SPBalanceWeights balanceWeights);
//		virtual void set_connection_weights(SmartPtr<ConnectionWeights<dim> >);

		virtual ConstSPProcessHierarchy current_process_hierarchy() const;
		virtual ConstSPProcessHierarchy next_process_hierarchy() const;

		virtual bool supports_balance_weights() const;
		virtual bool supports_connection_weights() const;
		virtual bool supports_repartitioning() const			{return false;}

		virtual bool partition(size_t baseLvl, size_t elementThreshold);

		virtual SubsetHandler& get_partitions();
		virtual const std::vector<int>* get_process_map() const;

	/**	static partitioning should be used if a globally refined multigrid is
	 * is used throughout a simulation. Static partitioning drastically reduces
	 * global communication. The drawback is, that static partitioning e.g.
	 * can't be used to properly rebalance a distributed adaptive grid.*/
		void enable_static_partitioning(bool enable);
		bool static_partitioning_enabled() const;

	private:
		enum constants{
			UNCLASSIFIED = 0,
			LEFT = 1,
			RIGHT = 1 << 1,
			CUTTING = LEFT | RIGHT,
			CUTTING_CENTER_LEFT,
			CUTTING_CENTER_RIGHT,
			TOTAL,
			NUM_CONSTANTS
		};

		static const size_t s_invalidIndex = -1;

		struct Entry{
			elem_t* elem;
			size_t	next;
			Entry(elem_t* e) : elem(e), next(s_invalidIndex)	{}
		};

		struct ElemList{
			ElemList() :
				m_entries(NULL), m_first(s_invalidIndex), m_last(s_invalidIndex), m_num(0)		{}
			ElemList(std::vector<Entry>* entries) :
				m_entries(entries), m_first(s_invalidIndex), m_last(s_invalidIndex), m_num(0)	{}

			void set_entry_list(std::vector<Entry>* entries)	{m_entries = entries;}

			void add(size_t entryInd)
			{
				if(m_first == s_invalidIndex){
					m_first = m_last = entryInd;
					m_num = 1;
				}
				else{
					(*m_entries)[m_last].next = entryInd;
					m_last = entryInd;
					++m_num;
				}
				(*m_entries)[entryInd].next = s_invalidIndex;
			}

			void clear()
			{
				m_first = m_last = s_invalidIndex;
				m_num = 0;
			}

			size_t size() const	{return m_num;}
			bool empty() const {return m_num == 0;}

			size_t first() const				{return m_first;}
			size_t last() const					{return m_last;}
			size_t next(size_t entryInd) const	{return (*m_entries)[entryInd].next;}
			elem_t* elem(size_t entryInd) const	{return (*m_entries)[entryInd].elem;}

			std::vector<Entry>* entries()		{return m_entries;}

			private:
				std::vector<Entry>* m_entries;
				size_t m_first;
				size_t m_last;
				size_t m_num;
		};


		struct TreeNode{
			ElemList	elems;
			int			firstProc;
			int			numTargetProcs;

			number		ratioLeft;
			size_t		firstChildNode;

			vector_t	center;
			vector_t	boxMin;
			vector_t	boxMax;
			number		totalWeight;

			number		splitValue;
			int			splitAxis;
			number		minSplitValue;
			number		maxSplitValue;

			bool		bisectionComplete;
		};


		void perform_bisection(int numTargetProcs, int minLvl, int maxLvl,
							   int partitionLvl, ANumber aWeight,
							   pcl::ProcessCommunicator com);

		void control_bisection(ISubsetHandler& partitionSH,
							   std::vector<TreeNode>& treeNodes, ANumber aWeight,
							   number maxChildWeight, pcl::ProcessCommunicator& com);

	/** if splitAxis is set to -1, then the split-axis is chosen for each node as
	 * the axis along which the geometry of the node has the longest extension.*/
		void bisect_elements(std::vector<TreeNode>& childNodesOut,
							 std::vector<TreeNode>& parentNodes,
							 ANumber aWeight, number maxChildWeight,
							 pcl::ProcessCommunicator& com,
							 int cutRecursion, int splitAxis);

	///	returns the next valid split axis to m_lastSplitAxis.
	/** The method starts from m_startSplitAxis, cycling through all valid split axis.
	 * \note	The method uses m_lastSplitAxis to determine the last used axis.
	 *			It also sets m_lastSplitAxis to the next axis.
	 *			In order to force the method to restart from m_startSplitAxis you may
	 *			set m_lastSplitAxis to -1.*/
		int get_next_split_axis();

	///	returns the next valid split axis to the specified one
		int get_next_split_axis(int lastAxis) const;


		void copy_partitions_to_children(ISubsetHandler& partitionSH, int lvl);

		void calculate_global_dimensions(std::vector<TreeNode>& treeNodes,
										 number maxChildWeight, ANumber aWeight,
										 pcl::ProcessCommunicator& com);

		void gather_weights_from_level(int baseLvl, int childLvl, ANumber aWeight,
									   bool copyToVMastersOnBaseLvl, bool markedElemsOnly);

		int classify_elem(elem_t* e, int splitAxis, number splitValue);

		void improve_split_values(std::vector<TreeNode>& treeNodes,
								  size_t maxIterations, ANumber aWeight,
								  pcl::ProcessCommunicator& com);
	

		MultiGrid*								m_mg;
		apos_t									m_aPos;
		aapos_t									m_aaPos;
		SmartPtr<SubsetHandler>					m_sh;
		SPProcessHierarchy						m_processHierarchy;
		SPProcessHierarchy						m_nextProcessHierarchy;
		pcl::InterfaceCommunicator<layout_t>	m_intfcCom;
		std::vector<Entry>						m_entries;

		SPBalanceWeights						m_balanceWeights;

		bool	m_staticPartitioning;

		number	m_tolerance;
		size_t	m_splitImproveIterations;

		int		m_longestSplitAxisEnabled;
		int		m_startSplitAxis;
		int		m_lastSplitAxis;
		bool	m_splitAxisEnabled[dim];
		int		m_numSplitAxisEnabled;
		int		m_firstSplitAxisEnabled;

	//	the following parameters are only used if the partitioner operates in
	//	static mode
		std::vector<int>	m_procMap;
		int 				m_highestRedistLevel;
};

///	\}

}// end of namespace

#endif
