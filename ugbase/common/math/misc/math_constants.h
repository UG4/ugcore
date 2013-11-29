/*
 * math_constants.h
 *
 *  Created on: 29.11.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UGMATH__MATH_CONSTANTS_H_
#define __H__UGMATH__MATH_CONSTANTS_H_

#include <cmath> //math.h>
#include "common/types.h"

namespace ug
{
	const number SMALL = 1.0e-12;
	const number SMALL_SQ = SMALL * SMALL;

	const number PI = M_PI;

}//	end of namespace

#endif /* __H__UGMATH__MATH_CONSTANTS_H_ */
