/*
 * profile_call.h
 *
 *  Created on: 23.07.2013
 *      Author: mrupp
 */



#ifndef PROFILE_CALL_H_
#define PROFILE_CALL_H_

#ifdef SHINY_CALL_LOGGING

#include <vector>

namespace ug{

struct ProfileCall
{
	ProfileCall(Shiny::ProfileNode *_p)
	{
		p = _p;
		c = clock();
	}
	Shiny::ProfileNode *p;
	clock_t c;
};

extern std::vector<ProfileCall> profileCalls;

}

#endif


#endif /* PROFILE_CALL_H_ */


