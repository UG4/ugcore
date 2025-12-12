/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Lisa Grau, Andreas Vogel
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


#ifndef __H__UG__LIB_DISC__QUADRATURE__NEWTON_COTES__
#define __H__UG__LIB_DISC__QUADRATURE__NEWTON_COTES__

#include "../quadrature.h"

namespace ug {


/**
 * This class provides Newton-Cotes integrals for the 1d line [0,1]. See e.g.
 * https://en.wikipedia.org/wiki/Newton%E2%80%93Cotes_formulas for details.
 *
 * The implemented rules are auto-generate using mathematica. If higher orders
 * are needed, rerun the corresponding file.
 */
class NewtonCotes : public QuadratureRule<1>
{
	public:
	/// constructor
		explicit NewtonCotes(size_t order);

	/// destructor
		~NewtonCotes() override;
};

} // namespace ug

#endif