
#ifndef __H__UG__DEBUG_PRINT_H_
#define __H__UG__DEBUG_PRINT_H_

#include <string>
#include <sstream>
#include "error.h"

namespace ug{

template<typename T>
void PrintVector(const T &v, std::string desc="")
{
	UG_LOG("================ " << desc << " ==================\n");
	for(size_t i=0; i<v.size(); i++)
	{
		UG_LOG(i << " = " << v[i] << "\n");
	}
	UG_LOG("==================================================\n");
}

#define PRINT_VECTOR(v, msg) {std::stringstream ss; ss << msg; PrintVector(v, ss.str()); }

#define UG_LINE_LOG() UG_LOG("\n-------LINE " << __FILE__ << ":" << __LINE__ << " (" << __PRETTY_FUNCTION__ << ") -----------\n");
}

#endif /* __H__UG__DEBUG_PRINT_H_ */
