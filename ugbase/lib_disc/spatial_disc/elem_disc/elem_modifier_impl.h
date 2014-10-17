/*
 * elem_modifier_impl.h
 *
 *  Created on: 04.10.2014
 *      Author: suze
 */

#ifndef ELEM_MODIFIER_IMPL_
#define ELEM_MODIFIER_IMPL_

namespace ug{

template <typename TDomain>
void IElemDiscModifier<TDomain>::
preprocess(LocalVector& u, LocalVector& d, LocalVector& tmpD, GridObject* elem,
			MathVector<dim> vCornerCoords[], LocalIndices& ind)
{
	UG_THROW("'preprocess_def' not implemented!\n");
}

template <typename TDomain>
void IElemDiscModifier<TDomain>::
preprocess(LocalVector& u, LocalMatrix& J, GridObject* elem,
			MathVector<dim> vCornerCoords[], LocalIndices& ind)
{
	UG_THROW("'preprocess_jac' not implemented!\n");
}


template <typename TDomain>
void IElemDiscModifier<TDomain>::
postprocess(const LocalVector& u, LocalVector& d, LocalIndices& ind)
{
	UG_THROW("'postprocess_def' not implemented!\n");
}

template <typename TDomain>
void IElemDiscModifier<TDomain>::
postprocess(const LocalVector& u, LocalMatrix& J, LocalIndices& ind)
{
	UG_THROW("'postprocess_jac' not implemented!\n");
}

} // end name space ug

#endif /* ELEM_MODIFIER_IMPL_ */
