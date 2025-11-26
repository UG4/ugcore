/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG_BRIDGE__ERROR__
#define __H__UG_BRIDGE__ERROR__

#include "common/common.h"
#include "common/error.h"

namespace ug{
namespace bridge{

/// \addtogroup registry
/// \{

struct UGRegistryError : UGError
{
	UGRegistryError(const std::string &name_,
	                const std::string &msg_,
	                const char* file = " -- no file -- ",
	                const unsigned long line = 0)
	:	UGError(msg_, file, line),
		name(name_), msg(msg_)
	{
		UG_ERR_LOG("#### Registry ERROR ("<<name_<<"):"<<msg_<<"\n");
	}

	UGRegistryError(const std::string &msg_,
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

// end group registry
/// \}

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


#endif