// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 25, 2013 (d,m,y)

#ifndef __H__UG__partitioner_parmetis__
#define __H__UG__partitioner_parmetis__

#include "lib_disc/domain.h"
#include "parallel_grid_layout.h"
#include "load_balancer.h"
#include "pcl/pcl_interface_communicator.h"
#include "util/parallel_dual_graph.h"

extern "C" {
	#include "metis.h"
	#include "parmetis.h"
}

namespace ug{

/*
// helper function for better conversion from number to idx_t
// i.e. for handling of special case, where d > numeric_limits<idx_t>::max()
inline idx_t cast_to_idx_t(const number& d)
{
	return	d >= std::numeric_limits<idx_t>::max() ?
			std::numeric_limits<idx_t>::max()
			: (idx_t) d;
}


// helper function for better addition of idx_t values
// i.e. for handling possible overflows when adding up very large values
inline void add_b_to_a(idx_t& a, const idx_t& b)
{
	a = b > std::numeric_limits<idx_t>::max() - a ?
		std::numeric_limits<idx_t>::max()
		: a+b;
}
*/

/// \addtogroup lib_grid_parallelization_distribution
///	\{

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
		virtual ~ICommunicationWeights() {};

		/**
		 * Get the weight of a specific connection.
		 * @param conn pointer to the connection in question
		 * @return weight for the side
		 */
		virtual number get_weight(GridObject* conn) = 0;

		/**
		 * Whether or not the given connection is to be assigned another weight.
		 * If true, this weight can be obtained by get_weight(),
		 * if false, no such weight must be requested.
		 * @param conn pointer to the connection in question
		 * @return true iff get_weight() can provide a new weight for the connection
		 */
		virtual bool reweigh(GridObject* conn) {return false;}
};



/// Implementation of ICommunicationCostWeights that can put specific weights on specific subsets
/**
 * This class enables the user to specify weights for specific subsets.
 * A fortiori, the user can choose infinite weights on a subset, if for
 * some reason this subset needs to be protected from being part of the
 * partition border.
 */
class SubsetCommunicationWeights : public ICommunicationWeights
{
	public:
		/// constructor
		SubsetCommunicationWeights(SmartPtr<IDomain<> > spDom) : m_sh(spDom->subset_handler()) {};

		/// destructor
		virtual ~SubsetCommunicationWeights() {};

		/// weight definition
		void set_weight_on_subset(number weight, int si)
		{
			if (m_weightMap.find(si) != m_weightMap.end())
			{
				UG_LOG("Warning: Mapping weights to subset " << si << "more than once "
					   "(in SubsetCommunicationCostWeights::set_weight_on_subset).");
			}
			m_weightMap[si] = weight;
		}

		/* As it turns out, ParMetis cannot handle this and will produce
		 * strange segfaults, probably due to uncaught overflows.
		/// infinite weight definition
		void set_infinite_weight_on_subset(int si)
		{
			number weight = std::numeric_limits<number>::has_infinity ?
							std::numeric_limits<number>::infinity()
							: std::numeric_limits<number>::max();
			set_weight_on_subset(weight, si);
		}
		*/

		/// getting the weights (inherited from ICommunicationCostWeights)
		virtual number get_weight(GridObject* conn)
		{
			if (!this->m_sh.valid())
				UG_THROW("Subset handler must be assigned to SubsetCommunicationCostWeights before it is used!");

			// get the subset for the connecting elem
			int si = m_sh->get_subset_index(conn);

			// check whether it is in one of the subsets with defined weight
			if (reweigh(conn))
				return m_weightMap[si];

			UG_THROW("Requested weight for element is not available. Check availability by calling"
					 "bool reweigh() before number get_weight().");
		}

		/// checking whether weight is available
		virtual bool reweigh(GridObject* conn)
		{
			// get the subset for the connecting elem
			int si = m_sh->get_subset_index(conn);

			return m_weightMap.find(si) != m_weightMap.end();
		}

	private:
		SmartPtr<MGSubsetHandler> m_sh;

		// storage for subset weight information
		std::map<int, number> m_weightMap;
};



/// Implementation of ICommunicationCostWeights that protects vertices of specific subsets
/**
 * This PartitionWeighting sets out to protect specific subsets from having
 * vertices in any process boundaries.
 * To that end, the get_weight() method will check  whether the subsets of the associated
 * vertices of any connection-type element belong to any of the subsets that are to be
 * protected; if so, a user-specified weight will be put onto that connection.
 * If more than one protected subset is involved the maximum weight will be used.
 */
class ProtectSubsetVerticesCommunicationWeights : public ICommunicationWeights
{
	public:
		typedef MultiGrid::traits<Vertex>::secure_container vertex_list_type;

		static const int NO_WEIGHT = 0;

		/// constructor
		ProtectSubsetVerticesCommunicationWeights(SmartPtr<IDomain<> > spDom)
		: m_sh(spDom->subset_handler()), m_weight(NO_WEIGHT) {};

		/// destructor
		virtual ~ProtectSubsetVerticesCommunicationWeights() {};

		/// weight definition
		void set_weight_on_subset(number weight, int si)
		{
			if (m_weightMap.find(si) != m_weightMap.end())
			{
				UG_LOG("Warning: Mapping weights to subset " << si << "more than once "
					   "(in SubsetCommunicationCostWeights::set_weight_on_subset).");
			}
			if (weight <= 0.0)
			{
				UG_THROW("Weights must be strictly positive!")
			}

			m_weightMap[si] = weight;
		}

		/// getting the weights (inherited from ICommunicationCostWeights)
		virtual number get_weight(GridObject* conn)
		{
			if (m_weight != NO_WEIGHT)
				return m_weight;

			UG_THROW("Requested weight for element is not available. Check availability by calling"
					 "bool reweigh() before number get_weight().");
		}

		/// checking whether weight is available
		virtual bool reweigh(GridObject* conn)
		{
			// reset local weight value
			m_weight = (number) NO_WEIGHT;
			bool protect = false;

			// check subset handler is available
			if (!this->m_sh.valid())
			{
				UG_THROW("Subset handler must be assigned to ProtectSubsetVerticesCommunicationCostWeights "
						 "before it is used!");
			}

			// get associated vertices
			vertex_list_type vert_list;
			m_sh->grid()->associated_elements(vert_list, conn);

			// loop vertices
			for (size_t i = 0; i < vert_list.size(); i++)
			{
				// get vertex subset index
				int si = m_sh->get_subset_index(vert_list[i]);

				// check whether it is to be protected
				std::map<int, number>::iterator it = m_weightMap.find(si);
				if (it != m_weightMap.end())
				{
					protect = true;

					// assign weight
					if (it->second > m_weight)
						m_weight = it->second;
				}
			}

			return protect;
		}

	private:
		// subset handler
		SmartPtr<MGSubsetHandler> m_sh;

		// storage for subset weight information
		std::map<int, number> m_weightMap;

		// temporary weight value (calculated by reweigh() and used by get_weight())
		number m_weight;
};



///	Parallel bisection partitioner
/**	The partitioner can be used inside a LoadBalancer or separately. It can
 * operate on serial and parallel multigrids. It is based on METIS and PARMETIS.
 * Please see ug's licensing page for more information on the licenses for METIS
 * and PARMETIS. Note that a special license is required to use PARMETIS for
 * commercial simulations and products.
 */
template <int dim>
class Partitioner_Parmetis : public IPartitioner{

	public:
		typedef IPartitioner base_class;

		using base_class::verbose;

		Partitioner_Parmetis();
		virtual ~Partitioner_Parmetis();

		virtual void set_grid(MultiGrid* mg, Attachment<MathVector<dim> > aPos);
		virtual void set_next_process_hierarchy(SPProcessHierarchy procHierarchy);
		virtual void set_balance_weights(SPBalanceWeights balanceWeights);
		virtual void set_communication_weights(SmartPtr<ICommunicationWeights> commWeights);

		virtual ConstSPProcessHierarchy current_process_hierarchy() const;
		virtual ConstSPProcessHierarchy next_process_hierarchy() const;

		virtual bool supports_balance_weights() const;
		virtual bool supports_communication_weights() const;
		virtual bool supports_repartitioning() const			{return true;}

		virtual bool partition(size_t baseLvl, size_t elementThreshold);

		virtual SubsetHandler& get_partitions();
		virtual const std::vector<int>* get_process_map() const;

		idx_t edge_cut_on_lvl(size_t lvl);

		void set_child_weight(int w);
		void set_sibling_weight(int w);

	///	Weights the cost of communication versus redistribution.
	/**	values in the range from 0.000001 to 1000000. A low value means that
	 * communication time is considered low compared to redistribution time while
	 * a high value means the contrary. Default is 1000.*/
		void set_itr_factor(float itr);

	/**
	 * @brief set factor of allowed imbalance
	 * This factor is set as allowed imbalance for each constraint.
	 * From METIS manual:
	 * "This is an array of size ncon that specifies the allowed load imbalance tolerance
	 * for each constraint. For the ith partition and jth constraint the allowed weight is
	 * the ubvec[j]*tpwgts[i*ncon+j] fraction of the jth’s constraint total weight.
	 * The load imbalances must be greater than 1.0."
	 *
	 * @param imb imbalance factor (default: 1.05)
	 */
		void set_allowed_imbalance_factor(float imb);

	private:
	///	fills m_aNumChildren with child-counts from levels baseLvl to topLvl.
	/**	Elements in topLvl have child-count 0.*/
		// dead code?
		template <typename TElem>
		void accumulate_child_counts(int baseLvl, int topLvl, AInt aInt,
									 pcl::InterfaceCommunicator<typename GridLayoutMap::Types<TElem>::Layout::LevelLayout>* p_intfcCom,
									 bool copyToVMastersOnBaseLvl = false);

	/**	writes the number of children in childLvl of each element in baseLvl into
	 *	the given int-attachment of the baseLvl elements.*/
		template <typename TElem>
		void gather_child_weights_from_level(int baseLvl, int childLvl, Attachment<idx_t> aWeight,
											bool copyToVMastersOnBaseLvl,
											pcl::InterfaceCommunicator<typename GridLayoutMap::Types<TElem>::Layout::LevelLayout>* p_intfcCom);

	///	Fills the array vwgtsOut with multi-constraint weights for elements in base-lvl
	/**	vwgtsOut will be resized to (#elemsInBaseLvl * (maxLvl-baseLvl+1)) resulting
	 * in a multi-constraint weights vector for elements in base-level. The i-th constraint
	 * of an element represents the weight of its child-elements in the i-th level.*/
		template <class TIter, typename TElem>
		void fill_multi_constraint_vertex_weights(std::vector<idx_t>& vwgtsOut,
											 int baseLvl, int maxLvl, Attachment<idx_t> aWeight,
											 bool fillVMastersOnBaseLvl,
											 TIter baseElemsBegin, TIter baseElemsEnd,
											 int numBaseElements,
											 pcl::InterfaceCommunicator<typename GridLayoutMap::Types<TElem>::Layout::LevelLayout>* p_intfcCom);
	/// the actual partitioning functionality
		template <typename TElem>
		bool partition_with_elem_type(size_t baseLvl, size_t elementThreshold);

	/**	make sure that m_numChildren contains a valid number for all elements
	 * on the given level!*/
		template <typename TElem>
		idx_t partition_level_metis(int baseLvl, int maxLvl, int numTargetProcs,
				  	  	  	  	   pcl::InterfaceCommunicator<typename GridLayoutMap::Types<TElem>::Layout::LevelLayout>* p_intfcCom);

	/**	make sure that m_numChildren contains a valid number for all elements
	 * on the given level!
	 * The given procCom should contain all processes which are potentially
	 * involved on the given level (before or after distribution).*/
		template <typename TElem>
		idx_t partition_level_parmetis(int baseLvl, int maxLvl, int numTargetProcs,
									  const pcl::ProcessCommunicator& procComAll,
									  ParallelDualGraph<TElem, idx_t>& pdg,
									  pcl::InterfaceCommunicator<typename GridLayoutMap::Types<TElem>::Layout::LevelLayout>* p_intfcCom);

	///	find dimension of highest-dimensional element in the associated grid
		void find_max_elem_dim();

		MultiGrid* m_mg;
		SubsetHandler m_sh;
		Attachment<idx_t> m_aWeightChildren;
		//Grid::AttachmentAccessor<elem_t, Attachment<idx_t> > m_aaWeightChildren;
		SPBalanceWeights	m_balanceWeights;
		SmartPtr<ICommunicationWeights>	m_communicationWeights;
		SPProcessHierarchy	m_processHierarchy;
		SPProcessHierarchy	m_nextProcessHierarchy;

		std::vector<idx_t> m_vEdgeCut;

		int	m_childWeight;
		int	m_siblingWeight;
		float m_comVsRedistRatio;
		float m_imbFactor;

		GridBaseObjectId m_elemType;
};

///	\}

} // end of namespace

#endif
