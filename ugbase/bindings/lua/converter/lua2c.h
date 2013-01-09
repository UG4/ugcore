/* 
 * File:   lua2c.h
 * Author: mrupp
 *
 * Created on 4. Dezember 2012, 17:02
 */

#ifndef LUA2C_H
#define	LUA2C_H

#include <stdio.h>
#include <dlfcn.h>
#include <string>
namespace ug{
namespace bridge {


class LUA2C
{

	typedef bool (*LUA2C_Function)(double *, double *) ;
	
	void* libHandle;
public:
	std::string name, pDyn;
	LUA2C_Function f;
	
	LUA2C() { f= NULL; name = "uninitialized"; pDyn = ""; libHandle = NULL;}
	int iIn, iOut;
	
	bool is_valid(int in, int out)
	{
		return f != NULL && in == iIn && out == iOut;
	}
	
	bool create(const char *functionName);
	
	/*double operator() (double ret[1], double [iDim])
	{
		return f(r, x);
	}*/
	
	virtual ~LUA2C();
};


}
}
#endif	/* LUA2C_H */

