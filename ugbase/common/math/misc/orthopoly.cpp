/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Author: Dmitry Logashenko
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
#include <cmath>

// own header
#include "orthopoly.h"

namespace ug
{

/** computes the (unscaled) Legendre polynomials
 *
 * The polynomials are \f$L_2\f$-orthogonal on \f$[-1, 1]\f$. They satisfy
 * the recursion \f$P_0 (x) = 1\f$, \f$P_1 (1) = x\f$,
 * \f$P_k (x) = ((2 k - 1) x P_{k-1} (x) - (k - 1) P_{k-2} (x)) / k\f$.
 * The \f$L_2\f$-norm of \f$P_k\f$ is \f$\sqrt {2 / (2 k + 1)}\f$.
 */
number LegendrePoly (size_t k, number x)
{
	if (k == 0) return 1;
	else if (k == 1) return x;
	
	return ((2 * k - 1) * x * LegendrePoly (k-1, x) - (k - 1) * LegendrePoly (k - 2, x)) / k;
}

/** computes the scaled Legendre polynomials
 *
 * The Legendre polynomials are \f$L_2\f$-orthogonal on \f$[-1, 1]\f$.
 * They satisfy the recursion \f$P_0 (x) = 1\f$, \f$P_1 (1) = x\f$,
 * \f$(n+1) P_k (x) = (2 k - 1) x P_{k-1} (x) - (k - 1) P_{n-2} (x)\f$.
 * The \f$L_2\f$-norm of \f$P_k\f$ is \f$\sqrt {2 / (2 k + 1)}\f$.
 * This function returns \f$\sqrt {(2 k + 1) / 2} P_k (x)\f$.
 */
number ScaledLegendrePoly (size_t k, number x)
{
	return sqrt (((number) (2 * k + 1)) / 2) * LegendrePoly (k, x);
}

/** computes the (unscaled) Chebyshev polynomials of the first kind
 *
 * The polynomials are orthogonal on \f$[-1, 1]\f$ w.r.t. the scalar product
 * \f$ \int_{-1}^1 \phi (x) \cdot \psi (x) \frac{1}{\sqrt {1 - x^2}} \, dx \f$
 * They satisfy the recursion \f$T_0 (x) = 1\f$, \f$T_1 (1) = x\f$,
 * \f$T_k (x) = 2 x P_{k-1} (x) - P_{k-2} (x)\f$.
 * The corresponding norm of \f$T_k\f$ is \f$\tfrac{\pi}{2}\f$ for \f$k > 0\f$
 * and \f$\pi\f$ for \f$k = 0\f$.
 */
number Chebyshev1Poly (size_t k, number x)
{
	if (k == 0) return 1;
	else if (k == 1) return x;
	
	return 2 * x * Chebyshev1Poly (k-1, x) - Chebyshev1Poly (k - 2, x);
}

/** computes the scaled Chebyshev polynomials of the first kind
 *
 * The polynomials are orthogonal on \f$[-1, 1]\f$ w.r.t. the scalar product
 * \f$ \int_{-1}^1 \phi (x) \cdot \psi (x) \frac{1}{\sqrt {1 - x^2}} \, dx \f$
 * They satisfy the recursion \f$T_0 (x) = 1\f$, \f$T_1 (1) = x\f$,
 * \f$T_k (x) = 2 x T_{k-1} (x) - T_{k-2} (x)\f$.
 * The corresponding norm of \f$T_k\f$ is \f$\tfrac{\pi}{2}\f$ for \f$k > 0\f$
 * and \f$\pi\f$ for \f$k = 0\f$.
 *
 * This function returns \f$T_k (x)\f$ divided by its norm.
 */
number ScaledChebyshev1Poly (size_t k, number x)
{
	if (k == 0) return Chebyshev1Poly (0, x) / M_PI;
	
	return Chebyshev1Poly (k, x) * 2 / M_PI;
}

/** computes the (unscaled) Chebyshev polynomials of the second kind
 *
 * The polynomials are orthogonal on \f$[-1, 1]\f$ w.r.t. the scalar product
 * \f$ \int_{-1}^1 \phi (x) \cdot \psi (x) \sqrt {1 - x^2} \, dx \f$
 * They satisfy the recursion \f$U_0 (x) = 1\f$, \f$U_1 (1) = 2 x\f$,
 * \f$U_k (x) = 2 x U_{k-1} (x) - U_{k-2} (x)\f$.
 * The corresponding norm of \f$T_k\f$ is \f$\tfrac{\pi}{2}\f$.
 */
number Chebyshev2Poly (size_t k, number x)
{
	if (k == 0) return 1;
	else if (k == 1) return 2 * x;
	
	return 2 * x * Chebyshev2Poly (k-1, x) - Chebyshev2Poly (k - 2, x);
}

/** computes the scaled Chebyshev polynomials of the second kind
 *
 * The polynomials are orthogonal on \f$[-1, 1]\f$ w.r.t. the scalar product
 * \f$ \int_{-1}^1 \phi (x) \cdot \psi (x) \sqrt {1 - x^2} \, dx \f$
 * They satisfy the recursion \f$U_0 (x) = 1\f$, \f$U_1 (1) = 2 x\f$,
 * \f$U_k (x) = 2 x U_{k-1} (x) - U_{k-2} (x)\f$.
 * The corresponding norm of \f$T_k\f$ is \f$\tfrac{\pi}{2}\f$.
 */
number ScaledChebyshev2Poly (size_t k, number x)
{
	return Chebyshev2Poly (k, x) * 2 / M_PI;
}

} // namespace ug

/* End of File */
