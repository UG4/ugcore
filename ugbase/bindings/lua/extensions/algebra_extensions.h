/*
 * algebra_extensions.h
 *
 *  Created on: 04.03.2011
 *      Author: kosta
 */

#ifndef ALGEBRA_EXTENSIONS_H_
#define ALGEBRA_EXTENSIONS_H_

#include "ug_bridge/ug_bridge.h"

namespace ug{
namespace bridge{

bool RegisterAlgebraExtensions(Registry& reg, const char* parentGroup);

}
}

#endif /* ALGEBRA_EXTENSIONS_H_ */
