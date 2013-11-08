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

template <int dim>
class Partitioner_DynamicBisection : public IPartitioner<dim>{
	public:
		typedef IPartitioner<dim> 			base_class;
		typedef typename base_class::elem_t	elem_t;
		typedef MathVector<dim>				vector_t;
		typedef Attachment<vector_t>		apos_t;
		typedef Grid::VertexAttachmentAccessor<apos_t>	aapos_t;
		typedef typename GridLayoutMap::Types<elem_t>::Layout::LevelLayout	layout_t;

		Partitioner_DynamicBisection();
		virtual ~Partitioner_DynamicBisection();

		virtual void set_grid(MultiGrid* mg, Attachment<MathVector<dim> > aPos);
		virtual void set_next_process_hierarchy(SPProcessHierarchy procHierarchy);
		virtual void set_balance_weights(SmartPtr<BalanceWeights<dim> >);
		virtual void set_connection_weights(SmartPtr<ConnectionWeights<dim> >);

		virtual ConstSPProcessHierarchy current_process_hierarchy() const;
		virtual ConstSPProcessHierarchy next_process_hierarchy() const;

		virtual bool supports_balance_weights() const;
		virtual bool supports_connection_weights() const;
		virtual bool supports_repartitioning() const			{return false;}

		virtual number estimate_distribution_quality(std::vector<number>* pLvlQualitiesOut = NULL);

		virtual void partition(size_t baseLvl, size_t elementThreshold);

		virtual SubsetHandler& get_partitions();
		virtual const std::vector<int>* get_process_map() const;

	private:
		enum constants{
			UNCLASSIFIED = 0,
			LEFT = 1,
			RIGHT = 1 << 1,
			CUTTING = LEFT | RIGHT,
			TOTAL,
			NUM_CONSTANTS
		};

		void copy_partitions_to_children(ISubsetHandler& partitionSH, int lvl);
		void perform_bisection(int minLvl, int maxLvl, int partitionLvl);
		void perform_bisection(int numTargetProcs, int minLvl, int maxLvl,
							   int partitionLvl, AInt aChildCount,
							   pcl::ProcessCommunicator com);
		void control_bisection(ISubsetHandler& partitionSH, const std::vector<elem_t*>& elems,
							 int maxNumChildren, int numTargetProcs, int firstProc,
							 AInt aChildCount, pcl::ProcessCommunicator& com);

		void bisect_elements(std::vector<elem_t*>& elemsLeftOut,
							std::vector<elem_t*>& elemsRightOut,
							const std::vector<elem_t*>& elems, number ratioLeft,
							AInt aChildCount, int maxNumChildren,
							pcl::ProcessCommunicator& com, int cutRecursion);

		void calculate_global_dimensions(vector_t& centerOut, vector_t& boxMinOut,
										 vector_t& boxMaxOut, const std::vector<elem_t*>& elems,
										 int maxNumChildren, AInt aChildCount,
										 pcl::ProcessCommunicator& com);
		void gather_child_counts_from_level(int baseLvl, int childLvl, AInt aInt,
											bool copyToVMastersOnBaseLvl);

		int classify_elem(elem_t* e, int splitDim, number splitValue);

		number find_split_value(const std::vector<elem_t*>& elems, int splitDim,
								number splitRatio, number initialGuess,
								number minValue, number maxValue,
								size_t maxIterations, AInt aChildCount,
								pcl::ProcessCommunicator& com);


		MultiGrid*								m_mg;
		apos_t									m_aPos;
		aapos_t									m_aaPos;
		SubsetHandler 							m_sh;
		SPProcessHierarchy						m_processHierarchy;
		SPProcessHierarchy						m_nextProcessHierarchy;
		pcl::InterfaceCommunicator<layout_t>	m_intfcCom;

};
}// end of namespace

#endif
