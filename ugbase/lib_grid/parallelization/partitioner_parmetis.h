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

/// \addtogroup lib_grid_parallelization_distribution
///	\{

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
//		virtual void set_connection_weights(SmartPtr<ConnectionWeights<dim> > conWeights);

		virtual ConstSPProcessHierarchy current_process_hierarchy() const;
		virtual ConstSPProcessHierarchy next_process_hierarchy() const;

		virtual bool supports_balance_weights() const;
		virtual bool supports_connection_weights() const;
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
		void accumulate_child_counts(int baseLvl, int topLvl, AInt aInt,
									 bool copyToVMastersOnBaseLvl = false);

	/**	writes the number of children in childLvl of each element in baseLvl into
	 *	the given int-attachment of the baseLvl elements.*/
		void gather_child_counts_from_level(int baseLvl, int childLvl, AInt aInt,
											bool copyToVMastersOnBaseLvl);

	///	Fills the array vwgtsOut with multi-constraint weights for elements in base-lvl
	/**	vwgtsOut will be resized to (#elemsInBaseLvl * (maxLvl-baseLvl+1)) resulting
	 * in a multi-constraint weights vector for elements in base-level. The i-th constraint
	 * of an element represents the weight of its child-elements in the i-th level.*/
		template <class TIter>
		void fill_multi_constraint_vertex_weights(std::vector<idx_t>& vwgtsOut,
											 int baseLvl, int maxLvl, AInt aInt,
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
		AInt m_aNumChildren;
		Grid::AttachmentAccessor<elem_t, AInt>	m_aaNumChildren;
		SPBalanceWeights	m_balanceWeights;
//		SPConnectionWeights	m_connectionWeights;
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
