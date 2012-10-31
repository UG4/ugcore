/*
 * error.h
 *
 *  Created on: 29.05.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG_BRIDGE__ERROR__
#define __H__UG_BRIDGE__ERROR__

#include "common/common.h"
#include "common/error.h"

namespace ug{
namespace bridge{

struct UGRegistryError : public UGError
{
	UGRegistryError(std::string name_,
	                                     std::string msg_,
	                     		         const char* file = " -- no file -- ",
	                     		         const unsigned long line = 0)
	:	UGError(msg_, file, line),
		name(name_), msg(msg_)
	{
		UG_ERR_LOG("#### Registry ERROR ("<<name_<<"):"<<msg_<<"\n");
	}

	UGRegistryError(std::string msg_,
	                     		         const char* file = " -- no file -- ",
	                     		         const unsigned long line = 0)
	:	UGError(msg_, file, line),
		name("-- no name --"), msg(msg_)
	{
		UG_ERR_LOG("#### Registry ERROR:"<<msg_<<"\n");
	}

	std::string name;
	std::string msg;
};


} // end namespace bridge
} // end namespace ug

#define UG_THROW_REGISTRY_ERROR(cls,msg) \
	{ug_throw_error(); std::stringstream ss; ss << msg; \
	throw(ug::bridge::UGRegistryError((cls),ss.str(),\
	                                                    __FILE__,__LINE__));}

#define UG_THROW_REGISTRY_MSG(msg) \
	{ug_throw_error(); std::stringstream ss; ss << msg; \
	throw(ug::bridge::UGRegistryError(ss.str(),\
	                                                    __FILE__,__LINE__));}


#endif /* __H__UG_BRIDGE__ERROR__ */
