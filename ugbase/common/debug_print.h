/*
 * debug_print.h
 *
 *  Created on: 10.05.2013
 *      Author: mrupp
 */

#ifndef __H__UG__DEBUG_PRINT_H_
#define __H__UG__DEBUG_PRINT_H_

#include <string>
#include <sstream>

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

}

#endif /* __H__UG__DEBUG_PRINT_H_ */
