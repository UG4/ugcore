/*
 * lib_algebra_bridge.cpp
 *
 *  Created on: 11.10.2010
 *      Author: andreasvogel
 */

// extern headers
#include <iostream>
#include <sstream>
#include <string>

// bridge
#include "bridge/bridge.h"
#include "lib_algebra/algebra_type.h"
#include "lib_algebra/cpu_algebra_types.h"

using namespace std;

namespace ug
{
namespace bridge
{

bool RegisterAMG(Registry& reg, string parentGroup);
bool RegisterEigensolver(Registry& reg, string parentGroup);
bool RegisterAlgebraCommon(Registry& reg, string parentGroup);
bool RegisterPreconditioner(Registry& reg, string parentGroup);
bool RegisterSolver(Registry& reg, string parentGroup);

bool RegisterLibAlgebra(Registry& reg, string parentGroup)
{
	RegisterAlgebraCommon(reg, parentGroup);
	RegisterPreconditioner(reg, parentGroup);
	RegisterSolver(reg, parentGroup);

	RegisterAMG(reg, parentGroup);
	RegisterEigensolver(reg, parentGroup);

	return true;
}



} // end namespace bridge
} // end namespace ug
