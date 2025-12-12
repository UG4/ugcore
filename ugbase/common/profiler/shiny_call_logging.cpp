/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "shiny_call_logging.h"

#include <vector>

#include "src/ShinyManager.h"
#include "profile_call.h"

namespace ug {

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
		profileCalls.push_back(ProfileCall(nullptr, tnow));
	else
		if( CheckEnoughTimePassedToNow(callsOnHold[N-1], tnow))
	{
		for(size_t i=0; i<N; i++)
			profileCalls.push_back(callsOnHold[i]);
		callsOnHold.clear();

		profileCalls.push_back(ProfileCall(nullptr, tnow));
	}
	else
	{
		callsOnHold.resize(N-1);
	}

}

}
