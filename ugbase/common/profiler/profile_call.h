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
#include "src/ShinyTools.h"

namespace ug{

struct ProfileCall
{
	ProfileCall()
	{
		p = NULL;
		c = 0;
	}
	ProfileCall(const ProfileCall &other)
	{
		p = other.p;
		c = other.c;
	}
	ProfileCall(Shiny::ProfileNode *_p)
	{
		p = _p;
		Shiny::GetTicks(&c);
	}
	Shiny::ProfileNode *p;
	Shiny::tick_t c;
};

extern std::vector<ProfileCall> profileCalls;

void FinishShinyCallLogging();
}

#endif


#endif /* PROFILE_CALL_H_ */


