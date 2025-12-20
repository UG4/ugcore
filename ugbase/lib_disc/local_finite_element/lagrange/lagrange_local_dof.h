/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__LAGRANGE_LOCAL_DOF__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__LAGRANGE_LOCAL_DOF__

#include "lib_disc/local_finite_element/local_dof_set.h"

namespace ug {

/// returns number of DoFs on element type for order p
size_t LagrangeNumDoFs(ReferenceObjectID_t elem, size_t p);

///	returns number of DoFs Subelement for an element type and order p
size_t LagrangeNumDoFOnSub(ReferenceObjectID_t elem,
                           ReferenceObjectID_t sub, size_t p);


///////////////////////////////////////////////////////////////////////////////
// Lagrange LocalDoFSets
///////////////////////////////////////////////////////////////////////////////

/// Lagrange DoF Set
template <typename TRefElem>
class LagrangeLDS : public LocalDoFSet
{
	public:
	///	constructor
		explicit LagrangeLDS(size_t order = 1);

	///	sets the order
		void set_order(size_t order);

	///	returns the type of reference element
		[[nodiscard]] ReferenceObjectID_t roid() const {return TRefElem::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		[[nodiscard]] size_t num_dof() const {return m_vLocalDoF.size();};

	///	returns the number of DoFs on a sub-geometric object type
		[[nodiscard]] size_t num_dof(ReferenceObjectID_t roid) const;

	///	returns the dof storage
		[[nodiscard]] const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}

	///	returns if the local dof position are exact
		[[nodiscard]] bool exact_position_available() const {return true;};

	protected:
		size_t p;		///< order
		std::vector<LocalDoF> m_vLocalDoF; 	///< association to geom obj
};

} //namespace ug

#endif