/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#include "../util_overloaded.h"
#include "ug.h"
#include "bridge/bridge.h"
#include "common/math/misc/orthopoly.h"
#include "registry/registry.h"

#include <cstdlib>
#include <string>

using namespace std;

namespace ug{

namespace bridge{

/// \ingroup misc_bridge
/// \{
template<typename TRegistry=Registry>
void RegisterBridge_OrthoPoly_(TRegistry& reg, string parentGroup)
{
	string grp(parentGroup);
	grp.append("/OrthoPoly");

	reg.add_function("LegendrePoly", &LegendrePoly, grp,
	                 "n#x", "", "Returns the value of the n-th Legendre polynomial at x");
	reg.add_function("SqNormOfLegendrePoly", &SqNormOfLegendrePoly, grp,
	                 "n", "", "Returns the scalar square of the n-th Legendre polynomial");
	reg.add_function("NormalizedLegendrePoly", &NormalizedLegendrePoly, grp,
	                 "n#x", "", "Returns the value of the normalized n-th Legendre polynomial at x");

	reg.add_function("Chebyshev1Poly", &Chebyshev1Poly, grp,
	                 "n#x", "", "Returns the value of the n-th Chebyshev polynomial of the first kind at x");
	reg.add_function("SqNormOfChebyshev1Poly", &SqNormOfChebyshev1Poly, grp,
	                 "n", "", "Returns the scalar square of the n-th Chebyshev polynomial of the first kind");
	reg.add_function("NormalizedChebyshev1Poly", &NormalizedChebyshev1Poly, grp,
	                 "n#x", "", "Returns the value of the normalized n-th Chebyshev polynomial of the first kind at x");

	reg.add_function("Chebyshev2Poly", &Chebyshev2Poly, grp,
	                 "n#x", "", "Returns the value of the n-th Chebyshev polynomial of the second kind at x");
	reg.add_function("SqNormOfChebyshev2Poly", &SqNormOfChebyshev2Poly, grp,
	                 "n", "", "Returns the scalar square of the n-th Chebyshev polynomial of the second kind");
	reg.add_function("NormalizedChebyshev2Poly", &NormalizedChebyshev2Poly, grp,
	                 "n#x", "", "Returns the value of the normalized n-th Chebyshev polynomial of the second kind at x");

}

// end group util_bridge
/// \}

}// end of namespace bridge

// Instantiate templates.
UG_REGISTRY_DEFINE(RegisterBridge_OrthoPoly);

}// end of namespace ug
