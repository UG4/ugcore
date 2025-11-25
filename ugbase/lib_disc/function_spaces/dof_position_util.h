/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
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
///�\{
template <int dim>
bool InnerDoFPosition(std::vector<MathVector<dim> >& vPos, ReferenceObjectID roid,
                      const std::vector<MathVector<dim> >& vCornerCoord, const LFEID& lfeID);


template <typename TDomain>
bool InnerDoFPosition(std::vector<MathVector<TDomain::dim> >& vPos,
                      GridObject* elem, const TDomain& domain, const LFEID& lfeID);
/// \}

/**
 * Returns the global DoF position on an element.
 */
///�\{
template <int dim>
bool DoFPosition(std::vector<MathVector<dim> >& vPos, ReferenceObjectID roid,
                 const std::vector<MathVector<dim> >& vCornerCoord, const LFEID& lfeID);

template <typename TDomain>
bool DoFPosition(std::vector<MathVector<TDomain::dim> >& vPos,
                 GridObject* elem, const TDomain& domain, const LFEID& lfeID);
/// \}

/**
 * returns the shape function values at given global positions on an element
 */
/// \{
template <int dim>
void ShapesAtGlobalPosition(std::vector<std::vector<number> >& vvShape,
                           const std::vector<MathVector<dim> >& vGlobPos,
                           ReferenceObjectID roid,
                           const std::vector<MathVector<dim> >& vCornerCoord,
                           const LFEID& lfeID);
template <typename TDomain>
void ShapesAtGlobalPosition(std::vector<std::vector<number> >& vvShape,
                           const std::vector<MathVector<TDomain::dim> >& vGlobPos,
                           GridObject* elem, const TDomain& domain, const LFEID& lfeID);
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
 * @param vPosPair			the array of positions (to be filled)
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
 * @param vPosPair			the array of positions (to be filled)
 */
template <typename TDomain>
void ExtractPositions(ConstSmartPtr<TDomain> domain,
                      ConstSmartPtr<DoFDistribution> dd,
                      size_t fct,
                      std::vector<std::pair<MathVector<TDomain::dim>, size_t> >& vPosPair);

/**
 * Checks that DoF Positions equal whether they are extracted by a macroelement
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

#endif