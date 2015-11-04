
#include <vector>
#include "src/ShinyManager.h"
#include "profile_call.h"

namespace ug{

std::vector<ProfileCall> profileCalls;

int g_ShinyCallLoggingMaxFreq = 1; // minimum frequency is 100 Hz = 10 ms per call


std::vector<ProfileCall> callsOnHold;

void ShinyCallLoggingStart()
{
	callsOnHold.push_back(ProfileCall(Shiny::ProfileManager::instance._curNode));
}

bool CheckEnoughTimePassedToNow(ProfileCall &pc, Shiny::tick_t tnow)
{
	static Shiny::tick_t freq = Shiny::GetTickFreq();  // in hz

	if( (tnow - pc.c)  > freq/g_ShinyCallLoggingMaxFreq)
		return true;
	else return false;
}

void FinishShinyCallLogging()
{
	for(size_t i=0; i<callsOnHold.size(); i++)
		profileCalls.push_back(callsOnHold[i]);
	callsOnHold.clear();
}

void ShinyCallLoggingEnd()
{
	size_t N = callsOnHold.size();

	Shiny::tick_t tnow;
	Shiny::GetTicks(&tnow);

	if(N == 0)
		profileCalls.push_back(ProfileCall(NULL, tnow));
	else
		if( CheckEnoughTimePassedToNow(callsOnHold[N-1], tnow))
	{
		for(size_t i=0; i<N; i++)
			profileCalls.push_back(callsOnHold[i]);
		callsOnHold.clear();

		profileCalls.push_back(ProfileCall(NULL, tnow));
	}
	else
	{
		callsOnHold.resize(N-1);
	}

}

}
