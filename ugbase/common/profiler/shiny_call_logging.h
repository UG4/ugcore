/*
 * shiny_call_logging.h
 *
 *  Created on: 08.07.2013
 *      Author: mrupp
 */

#ifndef SHINY_CALL_LOGGING_H_
#define SHINY_CALL_LOGGING_H_

#ifdef SHINY_CALL_LOGGING

namespace ug{
void ShinyCallLoggingStart();
void ShinyCallLoggingEnd();
extern int g_ShinyCallLoggingMaxFreq;
}

#define PROFILE_LOG_CALL_START() ug::ShinyCallLoggingStart();
#define PROFILE_LOG_CALL_END() ug::ShinyCallLoggingEnd();

#else
#define PROFILE_LOG_CALL_START()
#define PROFILE_LOG_CALL_END()

#endif



#endif /* SHINY_CALL_LOGGING_H_ */
