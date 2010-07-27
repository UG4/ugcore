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
		static const int world_dim = TWorldDim;
		static const int dim = TRefElem::dim;

	public:
		ReferenceMapping()
		{}

		void update(const MathVector<world_dim>* corners)
		{}

		bool local_to_global(	const MathVector<dim> loc_pos,
								MathVector<world_dim>& glob_pos) const
		{return false;}

		bool jacobian_transposed(	const MathVector<dim> loc_pos,
									MathMatrix<dim, world_dim>& JT) const
		{return false;}

		bool jacobian_transposed_inverse(	const MathVector<dim> loc_pos,
											MathMatrix<world_dim, dim>& JTInv) const
		{return false;}

		bool jacobian_det(const MathVector<dim> loc_pos, number& det) const
		{return false;}
};


} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_ELEMENT_MAPPING__ */
