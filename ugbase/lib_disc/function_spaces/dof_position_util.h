/*
 * dof_position_util.h
 *
 *  Created on: 17.08.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__DOF_POSITION_UTIL__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__DOF_POSITION_UTIL__

#include "common/common.h"
#include "common/math/ugmath.h"
#include "common/util/smart_pointer.h"

namespace ug {

/**
 * extracts the positions of the degrees of freedom and stores them into the
 * passed vector at the position of the algebraic index corresponding to the
 * degree of freedom.
 *
 * @param domain		the underlying domain
 * @param dd			the dof distribution
 * @param vPos			the array of positions (to be filled)
 */
template<typename TDomain, typename TDD>
void ExtractPositions(ConstSmartPtr<TDomain> domain, ConstSmartPtr<TDD> dd,
                      std::vector<MathVector<TDomain::dim> >& vPos);

/**
 * extracts the positions of the degrees of freedom and stores them into the
 * passed vector at the position of the algebraic index corresponding to the
 * degree of freedom.
 *
 * @param u				the underlying grid function
 * @param vPos			the array of positions (to be filled)
 */
template<typename TFunction>
void ExtractPositions(const TFunction &u,
                      std::vector<MathVector<TFunction::domain_type::dim> >& vPos)
{
	typedef typename TFunction::domain_type domain_type;
	typedef typename TFunction::dof_distribution_type dof_distribution_type;

	ExtractPositions<domain_type, dof_distribution_type>(u.domain(),
	                                                     u.dof_distribution(),
	                                                     vPos);
}

/**
 * extracts the positions of the degrees of freedom and stores them together with
 * the index in a std::pair into the passed vector at the position of the
 * algebraic index corresponding to the degree of freedom.
 *
 * @param domain		the underlying domain
 * @param dd			the dof distribution
 * @param vPos			the array of positions (to be filled)
 */
template <typename TDomain, typename TDD>
void ExtractPositions(ConstSmartPtr<TDomain> domain,
                      ConstSmartPtr<TDD> dd,
                      std::vector<std::pair<MathVector<TDomain::dim>, size_t> >& vPosPair);

} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__DOF_POSITION_UTIL__ */
