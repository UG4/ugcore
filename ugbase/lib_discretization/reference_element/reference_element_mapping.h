/*
 * reference_element_mapping.h
 *
 *  Created on: 13.05.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_ELEMENT_MAPPING__
#define __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_ELEMENT_MAPPING__

#include <cassert>
#include <iostream>
#include "common/common.h"
#include "common/math/ugmath.h"
#include "lib_grid/lg_base.h"

namespace ug{

template <typename TRefElem, int TWorldDim>
class ReferenceMapping
{
	public:
	///	world dimension (range space dimension)
		static const int worldDim = TWorldDim;

	///	reference dimension (domain space dimension)
		static const int dim = TRefElem::dim;

	public:
	///	Constructor
		ReferenceMapping() {}

	///	refresh mapping for new set of corners
		void update(const MathVector<worldDim>* corners) {}

	///	map local coordinate to global coordinate
		void local_to_global(	const MathVector<dim> locPos,
								MathVector<worldDim>& globPos) const {}

	///	returns transposed of jacobian
		void jacobian_transposed(	const MathVector<dim> locPos,
									MathMatrix<dim, worldDim>& JT) const {}

	///	returns transposed of the inverse of the jacobian
		void jacobian_transposed_inverse(	const MathVector<dim> locPos,
											MathMatrix<worldDim, dim>& JTInv) const {}

	///	returns the determinate of the jacobian
		number jacobian_det(const MathVector<dim> locPos) const {return 0.0;}
};


} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_ELEMENT_MAPPING__ */
