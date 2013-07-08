/*
 * shiny_call_logging.h
 *
 *  Created on: 08.07.2013
 *      Author: mrupp
 */

#ifndef SHINY_CALL_LOGGING_H_
#define SHINY_CALL_LOGGING_H_

#ifdef SHINY_CALL_LOGGING
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
#define PROFILE_LOG_CALL_START() profileCalls.push_back(ProfileCall(Shiny::ProfileManager::instance._curNode));
#define PROFILE_LOG_CALL_END() profileCalls.push_back(ProfileCall(NULL));

#else
#define PROFILE_LOG_CALL_START()
#define PROFILE_LOG_CALL_END()

#endif


#endif /* SHINY_CALL_LOGGING_H_ */
