/*
 * util.h
 *
 *  Created on: 26.09.2013
 *      Author: mrupp
 */

#ifndef UTIL_H_
#define UTIL_H_

#include "common/util/typename.h"

namespace ug{

template<typename TTo, typename TFrom>
TTo &DownCast(const TTo &pfrom, TFrom &p)
{
	TTo *t = dynamic_cast<TTo*>(&p);
	UG_ASSERT(t, "could not downcast " << TypeName(pfrom) << " to " << TypeName(p));
	return *t;
}

template<typename TTo, typename TFrom>
TTo &DownCast(const TTo &pfrom, const TFrom &p)
{
	const TTo *t = dynamic_cast<const TTo*>(&p);
	UG_ASSERT(t, "could not downcast " << TypeName(pfrom) << " to " << TypeName(p));
	return *t;
}

}
#endif /* UTIL_H_ */
