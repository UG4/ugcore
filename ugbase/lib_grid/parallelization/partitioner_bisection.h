// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Mar 1, 2013 (d,m,y)

#ifndef __H__UG__partitioner_bisection__
#define __H__UG__partitioner_bisection__

#include "load_balancer.h"

namespace ug{

template <int dim>
class Partitioner_Bisection : public IPartitioner<dim>{
	public:
		typedef IPartitioner<dim> base_class;
		typedef typename base_class::elem_t	elem_t;

		Partitioner_Bisection();
		virtual ~Partitioner_Bisection();

		virtual void set_grid(MultiGrid* mg, Attachment<MathVector<dim> > aPos);
		virtual void set_process_hierarchy(SPProcessHierarchy procHierarchy);
		virtual void set_balance_weights(SmartPtr<BalanceWeights<dim> >);
		virtual void set_connection_weights(SmartPtr<ConnectionWeights<dim> >);

		virtual bool supports_balance_weights() const;
		virtual bool supports_connection_weights() const;
		virtual bool supports_repartitioning() const			{return false;}

		virtual number estimate_distribution_quality(std::vector<number>* pLvlQualitiesOut = NULL);

		virtual void partition(size_t baseLvl, size_t elementThreshold);

		virtual SubsetHandler& get_partitions();
		virtual const std::vector<int>* get_process_map() const;

	private:
		MultiGrid* m_mg;
		Attachment<MathVector<dim> > m_aPos;
		SubsetHandler 		m_sh;
		SPProcessHierarchy	m_processHierarchy;
		std::vector<int>	m_procMap;
		int m_highestRedistLevel;
};

}// end of namespace

#endif
