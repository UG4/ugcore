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
#include "lib_disc/dof_manager/dof_distribution.h"

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
template<typename TDomain>
void ExtractPositions(ConstSmartPtr<TDomain> domain, ConstSmartPtr<DoFDistribution> dd,
                      std::vector<MathVector<TDomain::dim> >& vPos,
                      const std::vector<int>* pvMapGlobalToPatch = NULL);

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
	ExtractPositions(u.domain(),u.dof_distribution(), vPos);
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
template <typename TDomain>
void ExtractPositions(ConstSmartPtr<TDomain> domain,
                      ConstSmartPtr<DoFDistribution> dd,
                      std::vector<std::pair<MathVector<TDomain::dim>, size_t> >& vPosPair);

/**
 * extracts the positions of the degrees of freedom of a component
 * and stores them together with the index in a std::pair into the passed
 * vector at the position of the algebraic index corresponding to the degree of
 * freedom.
 *
 * @param domain		the underlying domain
 * @param dd			the dof distribution
 * @param fct			the component of the trial space
 * @param vPos			the array of positions (to be filled)
 */
template <typename TDomain>
void ExtractPositions(ConstSmartPtr<TDomain> domain,
                      ConstSmartPtr<DoFDistribution> dd,
                      const size_t fct,
                      std::vector<std::pair<MathVector<TDomain::dim>, size_t> >& vPosPair);

} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__DOF_POSITION_UTIL__ */
