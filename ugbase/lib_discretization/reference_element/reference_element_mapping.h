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

/**
 * This class describes the mapping from a reference element into the real
 * (physical) world. The mapping is initialized by the physical positions of
 * the vertices of the real world element. The order of those points must be
 * given as indicated by the corresponding reference element.
 *
 * Let \f$R\f$ be the reference element and \f$T\f$ be the element. Then, the
 * reference mapping is a mapping:
 * \f[
 * 	\phi:	R \mapsto T
 * \f]
 *
 * \tparam	TRefElem		reference element
 * \tparam	TWorldDim		world dimension
 */
template <typename TRefElem, int TWorldDim>
class ReferenceMapping
{
	public:
	///	world dimension (range space dimension)
		static const int worldDim = TWorldDim;

	///	reference dimension (domain space dimension)
		static const int dim = TRefElem::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static const bool isLinear = false;

	public:
	///	Default Constructor
		ReferenceMapping();

	///	Constructor setting the corners of the element
		ReferenceMapping(const MathVector<worldDim>* vCorner);

	///	refresh mapping for new set of corners
		void update(const MathVector<worldDim>* vCorner);

	///	map local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>& globPos,
		                     const MathVector<dim> locPos) const;

	///	returns transposed of jacobian
		void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
		                         const MathVector<dim> locPos) const;

	///	returns transposed of the inverse of the jacobian
		void jacobian_transposed_inverse(MathMatrix<worldDim, dim>& JTInv,
		                                 const MathVector<dim> locPos) const;

	///	returns the determinate of the jacobian
		number jacobian_det(const MathVector<dim> locPos) const;
};


} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_ELEMENT_MAPPING__ */
