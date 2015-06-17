/* 
 * \file	lua2c.cpp
 * \author	Martin Rupp
 *
 * Created on 4. Dezember 2012, 17:02
 */

#ifndef __LUACompiler_H__
#define	__LUACompiler_H__

#include <stdio.h>
#include <string>
#include "common/util/dynamic_library_util.h"
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

