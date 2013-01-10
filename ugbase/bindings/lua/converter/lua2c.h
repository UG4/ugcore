/* 
 * File:   lua2c.h
 * Author: mrupp
 *
 * Created on 4. Dezember 2012, 17:02
 */

#ifndef LUA2C_H
#define	LUA2C_H

#include <stdio.h>
#include <string>
#include "common/util/dynamic_library_util.h"

namespace ug{
namespace bridge {


class LUA2C
{
	
private:
	typedef bool (*LUA2C_Function)(double *, double *) ;
	
	DynLibHandle m_libHandle;
	std::string m_pDyn;;
	

public:
	std::string m_name;
	LUA2C_Function m_f;
	int m_iIn, m_iOut;
	
	LUA2C()
	{ 
		m_f= NULL; 
		m_name = "uninitialized"; 
		m_pDyn = ""; 
		m_libHandle = NULL;
	}	
	
	int num_in() const
	{
		return m_iIn;
	}
	
	int num_out() const
	{
		return m_iOut;
	}
	
	const string &name() const
	{
		return m_name;
	}
	
	bool is_valid() const
	{
		return m_f != NULL;
	}
	
	bool create(const char *functionName);
	
	inline bool call(double *ret, double *in) const
	{
		UG_ASSERT(m_f != NULL, "function " << m_name << " not valid");
		return m_f(ret, in);
	}
	
	virtual ~LUA2C();
};


}
}
#endif	/* LUA2C_H */

