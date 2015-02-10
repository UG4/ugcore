// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 25, 2013 (d,m,y)

#ifndef __H__UG__partitioner_parmetis__
#define __H__UG__partitioner_parmetis__

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
template <int dim>
class ICommunicationWeights
{
	public:
		typedef typename GeomObjBaseTypeByDim<dim>::base_obj_type elem_t;
		typedef typename elem_t::side side_t;

		virtual ~ICommunicationWeights(){};
		//virtual void refresh_weights(int baseLevel) = 0;

		/**
		 * Get the weight of a specific connection.
		 * @param conn pointer to the connection in question
		 * @return weight for the side
		 */
		virtual number get_weight(side_t* conn) = 0;

		/**
		 * Whether or not the given connection is to be assigned another weight.
		 * If true, this weight can be obtained by get_weight(),
		 * if false, no such weight must be obtained.
		 * @param conn pointer to the connection in question
		 * @return true iff get_weight() can provide a new weight for the connection
		 */
		virtual bool reweigh(side_t* conn) {return true;}
};



/// Implementation of ICommunicationCostWeights that can put specific weights on specific subsets
/**
 * This class enables the user to specify weights for specific subsets.
 * A fortiori, the user can choose infinite weights on a subset, if for
 * some reason this subset needs to be protected from being part of the
 * partition border.
 */
template <typename TDomain>
class SubsetCommunicationWeights : public ICommunicationWeights<TDomain::dim>
{
	public:
		typedef typename GeomObjBaseTypeByDim<TDomain::dim>::base_obj_type elem_t;
		typedef typename elem_t::side side_t;

		/// constructor
		SubsetCommunicationWeights(SmartPtr<TDomain> spDom) : m_sh(spDom->subset_handler()) {};

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
		virtual number get_weight(side_t* conn)
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
		virtual bool reweigh(side_t* conn)
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
template <typename TDomain>
class ProtectSubsetVerticesCommunicationWeights : public ICommunicationWeights<TDomain::dim>
{
	public:
		typedef typename GeomObjBaseTypeByDim<TDomain::dim>::base_obj_type elem_t;
		typedef typename elem_t::side side_t;
		typedef MultiGrid::traits<Vertex>::secure_container vertex_list_type;

		static const int NO_WEIGHT = 0;

		/// constructor
		ProtectSubsetVerticesCommunicationWeights(SmartPtr<TDomain> spDom)
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
		virtual number get_weight(side_t* conn)
		{
			if (m_weight != NO_WEIGHT)
				return m_weight;

			UG_THROW("Requested weight for element is not available. Check availability by calling"
					 "bool reweigh() before number get_weight().");
		}

		/// checking whether weight is available
		virtual bool reweigh(side_t* conn)
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
		typedef typename GeomObjBaseTypeByDim<dim>::base_obj_type	elem_t;
		typedef typename elem_t::side		side_t;

		using base_class::verbose;

		Partitioner_Parmetis();
		virtual ~Partitioner_Parmetis();

		virtual void set_grid(MultiGrid* mg, Attachment<MathVector<dim> > aPos);
		virtual void set_next_process_hierarchy(SPProcessHierarchy procHierarchy);
		virtual void set_balance_weights(SPBalanceWeights balanceWeights);
		virtual void set_communication_weights(SmartPtr<ICommunicationWeights<dim> > commWeights);

		virtual ConstSPProcessHierarchy current_process_hierarchy() const;
		virtual ConstSPProcessHierarchy next_process_hierarchy() const;

		virtual bool supports_balance_weights() const;
		virtual bool supports_communication_weights() const;
		virtual bool supports_repartitioning() const			{return true;}

		virtual bool partition(size_t baseLvl, size_t elementThreshold);

		virtual SubsetHandler& get_partitions();
		virtual const std::vector<int>* get_process_map() const;

		void set_child_weight(int w);
		void set_sibling_weight(int w);

	///	Weights the cost of communication versus redistribution.
	/**	values in the range from 0.000001 to 1000000. A low value means that
	 * communication time is considered low compared to redistribution time while
	 * a high value means the contrary. Default is 1000.*/
		void set_itr_factor(float itr);

	private:
	///	fills m_aNumChildren with child-counts from levels baseLvl to topLvl.
	/**	Elements in topLvl have child-count 0.*/
		// dead code?
		void accumulate_child_counts(int baseLvl, int topLvl, AInt aInt,
									 bool copyToVMastersOnBaseLvl = false);

	/**	writes the number of children in childLvl of each element in baseLvl into
	 *	the given int-attachment of the baseLvl elements.*/
		void gather_child_weights_from_level(int baseLvl, int childLvl, Attachment<idx_t> aWeight,
											bool copyToVMastersOnBaseLvl);

	///	Fills the array vwgtsOut with multi-constraint weights for elements in base-lvl
	/**	vwgtsOut will be resized to (#elemsInBaseLvl * (maxLvl-baseLvl+1)) resulting
	 * in a multi-constraint weights vector for elements in base-level. The i-th constraint
	 * of an element represents the weight of its child-elements in the i-th level.*/
		template <class TIter>
		void fill_multi_constraint_vertex_weights(std::vector<idx_t>& vwgtsOut,
											 int baseLvl, int maxLvl, Attachment<idx_t> aWeight,
											 bool fillVMastersOnBaseLvl,
											 TIter baseElemsBegin, TIter baseElemsEnd,
											 int numBaseElements);

	/**	make sure that m_numChildren contains a valid number for all elements
	 * on the given level!*/
		void partition_level_metis(int baseLvl, int maxLvl, int numTargetProcs);

	/**	make sure that m_numChildren contains a valid number for all elements
	 * on the given level!
	 * The given procCom should contain all processes which are potentially
	 * involved on the given level (before or after distribution).*/
		void partition_level_parmetis(int baseLvl, int maxLvl, int numTargetProcs,
									  const pcl::ProcessCommunicator& procComAll,
									  ParallelDualGraph<elem_t, idx_t>& pdg);

		typedef typename GridLayoutMap::Types<elem_t>::Layout::LevelLayout	layout_t;

		MultiGrid* m_mg;
		SubsetHandler m_sh;
		Attachment<idx_t> m_aWeightChildren;
		Grid::AttachmentAccessor<elem_t, Attachment<idx_t> >	m_aaWeightChildren;
		SPBalanceWeights	m_balanceWeights;
		SmartPtr<ICommunicationWeights<dim> >	m_communicationWeights;
		SPProcessHierarchy	m_processHierarchy;
		SPProcessHierarchy	m_nextProcessHierarchy;
		pcl::InterfaceCommunicator<layout_t>	m_intfcCom;

		int	m_childWeight;
		int	m_siblingWeight;
		float m_comVsRedistRatio;
};

///	\}

} // end of namespace

#endif
