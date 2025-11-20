/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_DISC__OPERATOR__LINEAR_OPERATOR__LEVEL_PRECONDITIONER_INTERFACE_H__
#define __H__LIB_DISC__OPERATOR__LINEAR_OPERATOR__LEVEL_PRECONDITIONER_INTERFACE_H__


#include "lib_algebra/operator/interface/preconditioner.h"
#include "lib_grid/tools/grid_level.h"


namespace ug {


/// A preconditioner for the multi-grid context which is aware of the grid level it operates on.
template <typename TAlgebra>
class ILevelPreconditioner : public IPreconditioner<TAlgebra>
{
	public:
		/// constructor
		ILevelPreconditioner()
		: IPreconditioner<TAlgebra>() {}

		/// constructor with grid level
		ILevelPreconditioner(const GridLevel& gl)
		: IPreconditioner<TAlgebra>(), m_gl(gl) {}

		/// clone constructor
		ILevelPreconditioner(const ILevelPreconditioner& parent)
		: IPreconditioner<TAlgebra>(parent), m_gl(parent.m_gl) {}

		~ILevelPreconditioner() override = default;

		/// set the grid level
		void set_grid_level(const GridLevel& gl)
		{
			if (gl != m_gl)
			{
				m_gl = gl;
				grid_level_has_changed();
			}
		}

		/// response to change in grid level
		virtual void grid_level_has_changed() {};

	protected:
		GridLevel m_gl;
};

} // end namespace ug


#endif