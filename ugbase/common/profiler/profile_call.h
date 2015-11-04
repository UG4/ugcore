


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
	ProfileCall(Shiny::ProfileNode *_p, Shiny::tick_t now)
	{
		p = _p;
		c = now;
	}
	Shiny::ProfileNode *p;
	Shiny::tick_t c;
};

extern std::vector<ProfileCall> profileCalls;

void FinishShinyCallLogging();
}

#endif


#endif /* PROFILE_CALL_H_ */


