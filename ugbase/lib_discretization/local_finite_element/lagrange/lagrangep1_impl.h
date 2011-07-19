/*
 * lagrangep1.h
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LAGRANGEP1__LAGRANGEP1_IMPL__
#define __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LAGRANGEP1__LAGRANGEP1_IMPL__

#include "lagrangep1.h"

namespace ug{

template <typename TRefElem, int TOrder>
void
LagrangeP1<TRefElem,TOrder>::
shapes(shape_type* sOut, const position_type& x) const
{
//	loop shapes
	for(size_t sh = 0; sh < num_sh(); ++sh)
		sOut[sh] = shape(sh, x);
}

template <typename TRefElem, int TOrder>
typename LagrangeP1<TRefElem,TOrder>::grad_type
LagrangeP1<TRefElem,TOrder>::
grad(size_t i, const position_type& x) const
{
//	tmp gradient
	grad_type tmpGrad;

//	evaluate
	evaluate_grad(i, x, tmpGrad);

//	return by value
	return tmpGrad;
}

template <typename TRefElem, int TOrder>
void
LagrangeP1<TRefElem,TOrder>::
grads(grad_type* gOut, const position_type& x) const
{
//	loop shapes
	for(size_t sh = 0; sh < num_sh(); ++sh)
		evaluate_grad(sh, x, gOut[sh]);
}


} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LAGRANGEP1__LAGRANGEP1_IMPL__ */
