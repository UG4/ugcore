/*
 * p1conform_impl.h
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#ifndef P1CONFORM_IMPL_H_
#define P1CONFORM_IMPL_H_

namespace ug{

template <typename TRefElem>
const LocalDoFPattern<TRefElem>&
P1conform<TRefElem>::
local_dof_pattern() const
{
	return m_ElementDoFPattern;
}

template <typename TRefElem>
P1conform<TRefElem>&
P1conform<TRefElem>::
inst()
{
	static P1conform<TRefElem> instance;
	return instance;
}


template <typename TRefElem>
P1conform<TRefElem>::
P1conform()
{
	m_ElementDoFPattern.set_num_dofs(ROID_VERTEX, 1);
};


}

#endif /* P1CONFORM_IMPL_H_ */
