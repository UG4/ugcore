/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef __LUACompiler_H__
#define	__LUACompiler_H__

#include <stdio.h>
#include <string>
#include "common/util/dynamic_library_util.h"
#include "bindings/lua/lua_function_handle.h"

namespace ug{

class VMAdd;

namespace bridge {



class LUACompiler
{
	
private:
	typedef int (*LUA2C_Function)(double *, const double *) ;
	
	DynLibHandle m_libHandle;
	std::string m_pDyn;
	VMAdd* vm;

public:
	std::string m_name;
	LUA2C_Function m_f;
	int m_iIn, m_iOut;
	bool bInitialized;
	bool bVM;
	LUACompiler()
	{ 
		m_f= NULL; 
		m_name = "uninitialized"; 
		m_pDyn = ""; 
		m_libHandle = NULL;
		bInitialized = false;
		bVM = false;
		vm = NULL;
	}
	
	int num_in() const
	{
		return m_iIn;
	}
	
	int num_out() const
	{
		return m_iOut;
	}
	
	const std::string &name() const
	{
		return m_name;
	}
	
	bool is_valid() const
	{
		return bInitialized;
	}
	
	bool create(const char *functionName, LuaFunctionHandle* pHandle = NULL);
	bool createVM(const char *functionName, LuaFunctionHandle* pHandle = NULL);
	bool createC(const char *functionName, LuaFunctionHandle* pHandle = NULL);
	
	bool call(double *ret, const double *in) const;
	virtual ~LUACompiler();
};


}
}
#endif	/* __LUACompiler_H__ */

