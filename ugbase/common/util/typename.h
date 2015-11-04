
#ifndef __H__UTIL__TYPENAME_H_
#define __H__UTIL__TYPENAME_H_
#include <string>
#include "demangle.h"

namespace ug{
template<typename T>
inline std::string TypeName(const T &t)
{
	return demangle(typeid(t).name());
}

template<typename T>
inline std::string TypeName()
{
	return demangle(typeid(T).name());
}

}
#endif /* __H__UTIL__TYPENAME_H_ */
