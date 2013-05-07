/*
 * contact_interface.h
 *
 *  Created on: 06.05.2013
 *      Author: raphaelprohl
 */

#ifndef CONTACT_INTERFACE_H_
#define CONTACT_INTERFACE_H_

// other ug4 modules
#include "common/common.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"

namespace ug{

template <typename TDomain, typename TGridFunction>
class IContactDisc
	: public IElemDisc<TDomain>
{
	private:
	///	Base class type
		typedef IElemDisc<TDomain> base_type;

	///	own type
		typedef IContactDisc<TDomain, TGridFunction> this_type;

	public:
	///	Domain type
		typedef typename base_type::domain_type domain_type;

	///	World dimension
		static const int dim = base_type::dim;

	///	Position type
		typedef typename base_type::position_type position_type;

	///	base element type of associated domain
		typedef typename domain_traits<TDomain::dim>::
				geometric_base_object TBaseElem;

	public:
		IContactDisc(){};

	/// Virtual destructor
		virtual ~IContactDisc() {}

		//virtual void set_elem_disc(SmartPtr<IElemDisc<TDomain> > elem) = 0;

		virtual void contactForces(TGridFunction& force, const TGridFunction& u,
				std::vector<MultiIndex<2> > vActiveSet) = 0;

		//virtual void contactForces_elem(LocalVector& locForce,
		//		TBaseElem* elem, const LocalVector& locU) = 0;
};

} //end namespace ug

#endif /* CONTACT_INTERFACE_H_ */
