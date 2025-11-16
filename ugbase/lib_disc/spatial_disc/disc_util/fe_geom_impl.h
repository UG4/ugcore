/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#include "fe_geom.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// FEGeometry
////////////////////////////////////////////////////////////////////////////////

template <	typename TElem,	int TWorldDim,
			typename TTrialSpace, typename TQuadratureRule>
FEGeometry<TElem,TWorldDim,TTrialSpace,TQuadratureRule>::
FEGeometry()
: m_rQuadRule(Provider<quad_rule_type>::get()),
  m_rTrialSpace(Provider<trial_space_type>::get())
{
	//	evaluate local shapes and gradients
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t sh = 0; sh < nsh; ++sh)
		{
			m_vvShape[ip][sh] = m_rTrialSpace.shape(sh, m_rQuadRule.point(ip));
			m_rTrialSpace.grad(m_vvGradLocal[ip][sh], sh, m_rQuadRule.point(ip));
		}
}

template <	typename TElem,	int TWorldDim,
			typename TTrialSpace, typename TQuadratureRule>
void
FEGeometry<TElem,TWorldDim,TTrialSpace,TQuadratureRule>::
update_local(ReferenceObjectID roid, const LFEID& lfeID, size_t orderQuad)
{
	if(roid != geometry_traits<TElem>::REFERENCE_OBJECT_ID)
		UG_THROW("FEGeometry::update: Geometry only for "
				<<geometry_traits<TElem>::REFERENCE_OBJECT_ID<<", but "
				<<roid<<" requested.");

	if(orderQuad > m_rQuadRule.order())
		UG_THROW("FEGeometry::update: Geometry only for order "
				<< m_rQuadRule.order()<<", but order "<<
				orderQuad<<" requested.");
}

template <	typename TElem,	int TWorldDim,
			typename TTrialSpace, typename TQuadratureRule>
void
FEGeometry<TElem,TWorldDim,TTrialSpace,TQuadratureRule>::
update(GridObject* elem, const MathVector<worldDim>* vCorner,
       const LFEID& lfeID, size_t orderQuad)
{
//	check
	UG_ASSERT(orderQuad <= m_rQuadRule.order(), "Wrong order requested.");

	UG_ASSERT(dynamic_cast<TElem*>(elem) != nullptr, "Wrong element type.");
	TElem* pElem = static_cast<TElem*>(elem);

//	check if element changed
	if(pElem == m_pElem) return;
	else m_pElem = pElem;

//	update the mapping for the new corners
	m_mapping.update(vCorner);

//	compute global integration points
	m_mapping.local_to_global(&m_vIPGlobal[0], local_ips(), nip);

//	evaluate global data
	m_mapping.jacobian_transposed_inverse(&m_vJTInv[0], &m_vDetJ[0],
	                                      local_ips(), nip);

// 	compute global gradients
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t sh = 0; sh < nsh; ++sh)
			MatVecMult(m_vvGradGlobal[ip][sh],
			           m_vJTInv[ip], m_vvGradLocal[ip][sh]);
}

} // end namespace ug
