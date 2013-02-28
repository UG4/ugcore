// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 27, 2013

#ifndef __H__UG__domain_load_balancer__
#define __H__UG__domain_load_balancer__

#include "lib_grid/parallelization/load_balancer.h"
#include "grid_function_serializer.h"

namespace ug{

///	A small wrapper for LoadBalancer which adds comfort methods to balance and distribute domains.
template <class TDomain, class TGridFct>
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
				GeomObjAttachmentSerializer<VertexBase, typename TDomain::position_attachment_type>::
								create(*dom->grid(), dom->position_attachment()));

			base_class::add_serializer(
				SubsetHandlerSerializer::create(*dom->subset_handler()));
		}

	///	grid functions added through this method are distributed along with the domain.
	/**	Make sure that added grid functions are consistent when distribution is performed!*/
		void add_serializer(SmartPtr<TGridFct> gridFct)
		{
			base_class::add_serializer(GridFunctionSerializer<TGridFct>::create(gridFct));
		}

	private:
		SmartPtr<TDomain>	m_dom;
};

}// end of namespace

#endif
