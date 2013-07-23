/*
 * shiny_call_logging.cpp
 *
 *  Created on: 23.07.2013
 *      Author: mrupp
 */

#include <vector>
#include "src/ShinyManager.h"
#include "profile_call.h"

namespace ug{

std::vector<ProfileCall> profileCalls;

void ShinyCallLoggingStart()
{
	profileCalls.push_back(ProfileCall(Shiny::ProfileManager::instance._curNode));
}

void ShinyCallLoggingEnd()
{
	profileCalls.push_back(ProfileCall(NULL));
}

}
