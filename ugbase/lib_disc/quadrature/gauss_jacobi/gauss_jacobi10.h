/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Lisa Grau, Andreas Vogel
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

#include "../quadrature.h"

#ifndef __H__UG__LIB_DISC__QUADRATURE__GAUSS_JACOBI10__
#define __H__UG__LIB_DISC__QUADRATURE__GAUSS_JACOBI10__

namespace ug{

/**
 * This class provides GaussJacobi integrals up to order 70
 * with alpha = 1 and beta = 0. For further information see e.g.
 *
 * Rathod, Venkatesh, Gauss Legendre - Gauss Jacobi Quadrature Rules over
 * a Tetrahedral Region, Int. J. Math Analysis, Vol. 5, 2011 (4), 189-198
 *
 * J. Villadsen and M.L. Michelsen, Solution of differential equation models by
 * polynomial approximation, Prentice Hall Inc, Englewood Cliffs,
 * New Jersey 07632 (1978)
 */
class GaussJacobi10 : public QuadratureRule<1>
{
	public:
	///	constructor
		GaussJacobi10(size_t order);

	///	destructor
		~GaussJacobi10();
};

} // namespace ug

#endif