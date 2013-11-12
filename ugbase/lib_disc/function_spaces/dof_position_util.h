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
 * Returns the global DoF position on an element.
 */
///Ê\{
template <int dim>
bool InnerDoFPosition(std::vector<MathVector<dim> >& vPos, const ReferenceObjectID roid,
                      const std::vector<MathVector<dim> >& vCornerCoord, const LFEID& lfeID);


template <typename TDomain>
bool InnerDoFPosition(std::vector<MathVector<TDomain::dim> >& vPos,
                      GeometricObject* elem, const TDomain& domain, const LFEID& lfeID);
/// \}

/**
 * Returns the global DoF position on an element.
 */
///Ê\{
template <int dim>
bool DoFPosition(std::vector<MathVector<dim> >& vPos, const ReferenceObjectID roid,
                 const std::vector<MathVector<dim> >& vCornerCoord, const LFEID& lfeID);

template <typename TDomain>
bool DoFPosition(std::vector<MathVector<TDomain::dim> >& vPos,
                 GeometricObject* elem, const TDomain& domain, const LFEID& lfeID);
/// \}

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

/**
 * Checks that DoF Positions equal wheather they are extracted by a macroelement
 * or a subelement.
 */
/// \{
template<typename TDomain>
bool CheckDoFPositions(ConstSmartPtr<TDomain> domain, ConstSmartPtr<DoFDistribution> dd);

template<typename TFunction>
bool CheckDoFPositions(const TFunction &u)
{
	return CheckDoFPositions(u.domain(),u.dof_distribution());
}

template <typename TDomain>
void ExtractAlgebraIndices(ConstSmartPtr<TDomain> domain,
                      ConstSmartPtr<DoFDistribution> dd,
                      std::vector<size_t> &fctIndex);
/// \}
} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__DOF_POSITION_UTIL__ */
