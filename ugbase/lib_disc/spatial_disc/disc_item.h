/*
 * disc_item.h
 *
 *  Created on: 31.08.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DISC_ITEM__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DISC_ITEM__

namespace ug{

template <typename TDomain, typename TAlgebra>
class IDomainConstraint;

template <typename TDomain, typename TAlgebra>
class IDiscretizationItem
{
	public:
	///	Type of Domain
		typedef TDomain domain_type;

	///	Type of algebra
		typedef TAlgebra algebra_type;

	public:
	///	returns the number of element discs
		virtual size_t num_elem_disc() const = 0;

	///	returns the element disc
		virtual SmartPtr<IElemDisc<TDomain> > elem_disc(size_t i) = 0;

	///	returns the number of constraints
		virtual size_t num_constraint() const = 0;

	///	returns an element disc
		virtual SmartPtr<IDomainConstraint<TDomain, TAlgebra> > constraint(size_t i) = 0;

	///	virtual destructor
		virtual ~IDiscretizationItem() {}
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DISC_ITEM__ */
