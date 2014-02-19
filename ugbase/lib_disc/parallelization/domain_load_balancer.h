// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 27, 2013

#ifndef __H__UG__domain_load_balancer__
#define __H__UG__domain_load_balancer__

#include "lib_grid/parallelization/load_balancer.h"
#include "lib_grid/parallelization/load_balancer_util.h"

namespace ug{

///	Creates a process-hierarchy that fullfills the given conditions.
template <class TDomain>
SPProcessHierarchy
CreateProcessHierarchy(TDomain& dom, size_t minNumElemsPerProcPerLvl,
					   size_t maxNumRedistProcs, size_t maxNumProcs,
					   int minDistLvl, int maxLevelsWithoutRedist)
{
	const DomainInfo& domInf = dom.domain_info();
	std::vector<size_t> numElemsOnLvl;
	numElemsOnLvl.reserve(domInf.num_levels());
	for(size_t i = 0; i < domInf.num_levels(); ++i)
		numElemsOnLvl.push_back(domInf.num_elements_on_level(i));

	if(numElemsOnLvl.empty()){
		return ProcessHierarchy::create();
	}

	return CreateProcessHierarchy(&numElemsOnLvl.front(), numElemsOnLvl.size(),
								  minNumElemsPerProcPerLvl, maxNumRedistProcs,
								  maxNumProcs, minDistLvl, maxLevelsWithoutRedist);
}


///	A small wrapper for LoadBalancer which adds comfort methods to balance and distribute domains.
template <class TDomain>
class DomainLoadBalancer : public LoadBalancer<TDomain::dim>
{
	typedef LoadBalancer<TDomain::dim> base_class;

	public:
		using base_class::add_serializer;

	/**	Adds serializers for the domains position attachment and for its subset handler.*/
		DomainLoadBalancer(SmartPtr<TDomain> dom) : m_dom(dom)
		{
			base_class::set_grid(dom->grid().get(), dom->position_attachment());

			base_class::add_serializer(
				GeomObjAttachmentSerializer<Vertex, typename TDomain::position_attachment_type>::
								create(*dom->grid(), dom->position_attachment()));

			base_class::add_serializer(
				SubsetHandlerSerializer::create(*dom->subset_handler()));

			std::vector<std::string> additionalSHNames = dom->additional_subset_handler_names();
			for(size_t i = 0; i < additionalSHNames.size(); ++i){
				SmartPtr<ISubsetHandler> sh = dom->additional_subset_handler(additionalSHNames[i]);
				if(sh.valid()){
					base_class::add_serializer(SubsetHandlerSerializer::create(*sh));
				}
			}
		}

	private:
		SmartPtr<TDomain>	m_dom;
};

}// end of namespace

#endif
