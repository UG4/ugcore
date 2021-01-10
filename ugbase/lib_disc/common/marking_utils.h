/*
 * Copyright (c) 2011-2018:  G-CSC, Goethe University Frankfurt
 * Author: Markus Breit
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

#ifndef UGBASE__LIB_GRID__REFINEMENT__MARKING_UTILS_H
#define UGBASE__LIB_GRID__REFINEMENT__MARKING_UTILS_H


#include "common/types.h"               // for number
#include "common/util/smart_pointer.h"  // for SmartPtr, ConstSmartPtr

#include <string>
#include <vector>


namespace ug {


// forward declarations
class IRefiner;
template <typename TDomain> class ApproximationSpace;




/// Mark all surface elements for refinement.
template <typename TDomain>
void MarkGlobal(SmartPtr<IRefiner> refiner, SmartPtr<TDomain> domain);


/// Mark all surface elements of specific subsets for refinement.
template <typename TDomain>
void MarkSubsets
(
	SmartPtr<IRefiner> refiner,
	SmartPtr<TDomain> domain,
	const std::vector<std::string>& vSubset
);


/**
 * @brief Mark all anisotropic elements of the surface grid for ansiotropic refinement
 *
 * Whether an element is isotropic or not is decided using the is_anisotropic() functions
 * from lib_grid's anisotropy_util and the given threshold ratio.
 */
template <typename TDomain>
void MarkAnisotropic
(
	SmartPtr<IRefiner> refiner,
	SmartPtr<TDomain> domain,
	number thresholdRatio
);


/**
 * @brief Mark all elements that are anisotropic in direction of the x-axis for ansiotropic refinement
 *
 * Whether an element is isotropic or not is decided using the is_anisotropic() functions
 * from lib_grid's anisotropy_util and the given threshold ratio.
 *
 * "Anisotropic in direction of the x-axis" means that the long edges (only the first one is checked)
 * point more or less in x-direction, to be precise: the normalized vector connecting this edge's
 * vertices has an x-entry of more than 0.9.
 */
template <typename TDomain>
void MarkAnisotropicOnlyX
(
	SmartPtr<IRefiner> refiner,
	SmartPtr<TDomain> domain,
	number thresholdRatio
);



/**
 *	Marks for refinement all (full-dim) elements neighboring grid elements
 *	that contain a degree of freedom (Lagrangian) whose value is outside a
 *	given range.
 *
 *	This function can be used to adaptively refine geometries with Q1 shape
 *	functions, which can become negative, e.g., in diffusion problems.
 *
 * @param refiner     refiner for hanging node refinement
 * @param u           solution grid function
 * @param cmp         component to check
 * @param lowerBnd    lower bound
 * @param upperBnd    upper bound
 */
template <typename TGridFunction>
void MarkOutOfRangeElems
(
	SmartPtr<IRefiner> refiner,
	ConstSmartPtr<TGridFunction> u,
	size_t cmp,
	number lowerBnd,
	number upperBnd
);



} // namespace ug


#include "marking_utils_impl.h"


#endif // UGBASE__LIB_GRID__REFINEMENT__MARKING_UTILS_H
