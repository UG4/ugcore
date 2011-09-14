/*
 * disc_item.h
 *
 *  Created on: 31.08.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG_LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_ITEM__
#define __H__UG_LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_ITEM__

namespace ug{

template <typename TDoFDistribution, typename TAlgebra>
class IConstraint;

template <	typename TDomain,
			typename TDoFDistribution,
			typename TAlgebra>
class IDiscretizationItem
{
	public:
	///	Type of Domain
		typedef TDomain domain_type;

	///	Type of DoF Distribution
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;

	///	Type of algebra
		typedef TAlgebra algebra_type;

	public:
	///	returns the number of element discs
		virtual size_t num_elem_disc() const = 0;

	///	returns the element disc
		virtual IDomainElemDisc<TDomain>* get_elem_disc(size_t i) = 0;

	///	returns the number of constraints
		virtual size_t num_constraint() const = 0;

	///	returns an element disc
		virtual IConstraint<TDoFDistribution, TAlgebra>* get_constraint(size_t i) = 0;

	///	virtual destructor
		virtual ~IDiscretizationItem() {}
};

} // end namespace ug

#endif /* __H__UG_LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_ITEM__ */
