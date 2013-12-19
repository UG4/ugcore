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


std::vector<ProfileCall> callsOnHold;

void ShinyCallLoggingStart()
{
	callsOnHold.push_back(ProfileCall(Shiny::ProfileManager::instance._curNode));
}

bool CheckEnoughTimePassedToNow(ProfileCall &pc)
{
	Shiny::tick_t tnow;
	Shiny::GetTicks(&tnow);
	static Shiny::tick_t freq = Shiny::GetTickFreq();  // in hz
	const size_t maxFreq = 10000; // minimum frequency is 10000 Hz = 0.1 ms per call
	if( (tnow - pc.c)  > freq/maxFreq)
		return true;
	else return false;
}

void FinishShinyCallLogging()
{
	for(size_t i=0; i<callsOnHold.size(); i++)
		profileCalls.push_back(callsOnHold[i]);
}

void ShinyCallLoggingEnd()
{
	size_t N = callsOnHold.size();

	if( CheckEnoughTimePassedToNow(callsOnHold[N-1]))
	{
		for(size_t i=0; i<N; i++)
			profileCalls.push_back(callsOnHold[i]);
		callsOnHold.clear();

		profileCalls.push_back(ProfileCall(NULL));
	}
	else
	{
		callsOnHold.resize(N-1);
	}

}

}
