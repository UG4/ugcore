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

#include "orthopoly.h"

#include <cmath>
#include "math_constants.h"



namespace ug {

/** returns the values of the Legendre polynomials
 *
 * The polynomials are \f$L_2\f$-orthogonal on \f$[-1, 1]\f$. They satisfy
 * the recursion \f$P_0 (x) = 1\f$, \f$P_1 (1) = x\f$,
 * \f$P_k (x) = ((2 k - 1) x P_{k-1} (x) - (k - 1) P_{k-2} (x)) / k\f$.
 * The \f$L_2\f$-norm of \f$P_k\f$ is \f$\sqrt {2 / (2 k + 1)}\f$.
 */
number LegendrePoly
(
	size_t k, ///< index of the polynomial, \f$k \ge 0\f$
	number x ///< argument of the polynomial
)
{
	if (k == 0) return 1;
	else if (k == 1) return x;
	
	return ((2 * k - 1) * x * LegendrePoly (k-1, x) - (k - 1) * LegendrePoly (k - 2, x)) / k;
}

/** returns the scalar square of the Legendre polynomials (the squared weighted norm)
 *
 * The polynomials are \f$L_2\f$-orthogonal on \f$[-1, 1]\f$. They satisfy
 * the recursion \f$P_0 (x) = 1\f$, \f$P_1 (1) = x\f$,
 * \f$P_k (x) = ((2 k - 1) x P_{k-1} (x) - (k - 1) P_{k-2} (x)) / k\f$.
 * The \f$L_2\f$-square of \f$P_k\f$ is \f$2 / (2 k + 1)\f$.
 * Note that this function returns not the \f$L_2\f$ norm but the weighted
 * norm \f$\sqrt {\frac{1}{b-a} \int_a^b P_k^2 (x) \, dx}\f$, where
 * \f$a=-1\f$, \f$b=1\f$.
 */
number SqNormOfLegendrePoly
(
	size_t k ///< index of the polynomial, \f$k \ge 0\f$
)
{
	return (number) 1 / (2 * k + 1);
}

/** returns the values of the normalized Legendre polynomials
 *
 * The Legendre polynomials are \f$L_2\f$-orthogonal on \f$[-1, 1]\f$.
 * They satisfy the recursion \f$P_0 (x) = 1\f$, \f$P_1 (1) = x\f$,
 * \f$(n+1) P_k (x) = (2 k - 1) x P_{k-1} (x) - (k - 1) P_{n-2} (x)\f$.
 * The \f$L_2\f$-norm of \f$P_k\f$ is \f$\sqrt {2 / (2 k + 1)}\f$.
 * This function returns \f$\sqrt {(2 k + 1) / 2} P_k (x)\f$.
 */
number NormalizedLegendrePoly
(
	size_t k, ///< index of the polynomial, \f$k \ge 0\f$
	number x ///< argument of the polynomial
)
{
	return sqrt (((number) (2 * k + 1))) * LegendrePoly (k, x);
}

/** returns the values of the Chebyshev polynomials of the first kind
 *
 * The polynomials are orthogonal on \f$[-1, 1]\f$ w.r.t. the scalar product
 * \f$ \int_{-1}^1 \phi (x) \cdot \psi (x) \frac{1}{\sqrt {1 - x^2}} \, dx \f$
 * They satisfy the recursion \f$T_0 (x) = 1\f$, \f$T_1 (1) = x\f$,
 * \f$T_k (x) = 2 x P_{k-1} (x) - P_{k-2} (x)\f$.
 * The corresponding norm of \f$T_k\f$ is \f$\tfrac{\pi}{2}\f$ for \f$k > 0\f$
 * and \f$\pi\f$ for \f$k = 0\f$.
 */
number Chebyshev1Poly
(
	size_t k, ///< index of the polynomial, \f$k \ge 0\f$
	number x ///< argument of the polynomial
)
{
	if (k == 0) return 1;
	else if (k == 1) return x;
	
	return 2 * x * Chebyshev1Poly (k-1, x) - Chebyshev1Poly (k - 2, x);
}

/** returns the scalar square of the Chebyshev polynomials of the first kind (the squared norm)
 *
 * The polynomials are orthogonal on \f$[-1, 1]\f$ w.r.t. the scalar product
 * \f$ \int_{-1}^1 \phi (x) \cdot \psi (x) \frac{1}{\sqrt {1 - x^2}} \, dx \f$
 * They satisfy the recursion \f$T_0 (x) = 1\f$, \f$T_1 (1) = x\f$,
 * \f$T_k (x) = 2 x P_{k-1} (x) - P_{k-2} (x)\f$.
 * The corresponding norm of \f$T_k\f$ is \f$\tfrac{\pi}{2}\f$ for \f$k > 0\f$
 * and \f$\pi\f$ for \f$k = 0\f$.
 */
number SqNormOfChebyshev1Poly
(
	size_t k ///< index of the polynomial, \f$k \ge 0\f$
)
{
	if (k == 0) return PI * PI;
	return PI * PI / 4;
}

/** returns the values the normalized Chebyshev polynomials of the first kind
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
number NormalizedChebyshev1Poly
(
	size_t k, ///< index of the polynomial, \f$k \ge 0\f$
	number x ///< argument of the polynomial
)
{
	if (k == 0) return Chebyshev1Poly (0, x) / PI;
	
	return Chebyshev1Poly (k, x) * 2 / PI;
}

/** returns the values of the Chebyshev polynomials of the second kind
 *
 * The polynomials are orthogonal on \f$[-1, 1]\f$ w.r.t. the scalar product
 * \f$ \int_{-1}^1 \phi (x) \cdot \psi (x) \sqrt {1 - x^2} \, dx \f$
 * They satisfy the recursion \f$U_0 (x) = 1\f$, \f$U_1 (1) = 2 x\f$,
 * \f$U_k (x) = 2 x U_{k-1} (x) - U_{k-2} (x)\f$.
 * The corresponding norm of \f$T_k\f$ is \f$\tfrac{\pi}{2}\f$.
 */
number Chebyshev2Poly
(
	size_t k, ///< index of the polynomial, \f$k \ge 0\f$
	number x ///< argument of the polynomial
)
{
	if (k == 0) return 1;
	else if (k == 1) return 2 * x;
	
	return 2 * x * Chebyshev2Poly (k-1, x) - Chebyshev2Poly (k - 2, x);
}

/** returns the scalar square of the Chebyshev polynomials of the second kind (the squared norm)
 *
 * The polynomials are orthogonal on \f$[-1, 1]\f$ w.r.t. the scalar product
 * \f$ \int_{-1}^1 \phi (x) \cdot \psi (x) \sqrt {1 - x^2} \, dx \f$
 * They satisfy the recursion \f$U_0 (x) = 1\f$, \f$U_1 (1) = 2 x\f$,
 * \f$U_k (x) = 2 x U_{k-1} (x) - U_{k-2} (x)\f$.
 * The corresponding norm of \f$T_k\f$ is \f$\tfrac{\pi}{2}\f$.
 */
number SqNormOfChebyshev2Poly
(
	size_t k ///< index of the polynomial, \f$k \ge 0\f$
)
{
	return PI * PI / 4;
}

/** returns the values the normalized Chebyshev polynomials of the second kind
 *
 * The polynomials are orthogonal on \f$[-1, 1]\f$ w.r.t. the scalar product
 * \f$ \int_{-1}^1 \phi (x) \cdot \psi (x) \sqrt {1 - x^2} \, dx \f$
 * They satisfy the recursion \f$U_0 (x) = 1\f$, \f$U_1 (1) = 2 x\f$,
 * \f$U_k (x) = 2 x U_{k-1} (x) - U_{k-2} (x)\f$.
 * The corresponding norm of \f$T_k\f$ is \f$\tfrac{\pi}{2}\f$.
 */
number NormalizedChebyshev2Poly
(
	size_t k, ///< index of the polynomial, \f$k \ge 0\f$
	number x ///< argument of the polynomial
)
{
	return Chebyshev2Poly (k, x) * 2 / PI;
}

} // namespace ug

