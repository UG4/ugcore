/*
 * ostream_util.h
 *
 *  Created on: 26.04.2013
 *      Author: mrupp
 */

#ifndef __H__UG_OSTREAM_UTIL_H_
#define __H__UG_OSTREAM_UTIL_H_

namespace ug
{

/// \addtogroup ugbase_common_io
/// \{

///used to unset flags set by std::scientific or std::fixed.
/**
 * example:
 * std::cout << std::scienfitic << myDouble << "\n";
 * std::cout << reset_floats << myPrecentage << "\n";
 * also in UG_LOG:
 * UG_LOG(reset_floats << myValue);
 * \note use ug::reset_floats if not in namespace ug.
 */
static inline std::ios_base&
reset_floats(std::ios_base& o)
{
	o.unsetf(std::ios_base::floatfield);
	return o;
}

// end group ugbase_common_io
/// \}

}
#endif /* __H__UG_OSTREAM_UTIL_H_ */
