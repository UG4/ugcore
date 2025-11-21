/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#include "profiler.h"

#ifdef UG_PROFILER

void ProfileNodeManager::
add(AutoProfileNode* node)
{
	inst().m_nodes.push(node);
}

void ProfileNodeManager::
release_latest()
{
	if(!inst().m_nodes.empty()){
		AutoProfileNode* node = inst().m_nodes.top();
		inst().m_nodes.pop();
		node->release();
	}
}


ProfileNodeManager::
~ProfileNodeManager()
{
//	release and deactivate all nodes
	while(!inst().m_nodes.empty())
		inst().release_latest();
}

ProfileNodeManager& ProfileNodeManager::
inst()
{
	static ProfileNodeManager pnm;
	return pnm;
}




#ifdef UG_PROFILER_SHINY
AutoProfileNode::AutoProfileNode() : m_bActive(true)
#endif
#if defined(UG_PROFILER_SCALASCA) || defined(UG_PROFILER_VAMPIR)
AutoProfileNode::AutoProfileNode(const char* name) : m_bActive(true), m_pName(name)
#endif
#ifdef UG_PROFILER_SCOREP
AutoProfileNode::AutoProfileNode(SCOREP_User_RegionHandle handle) : m_bActive(true), m_pHandle(handle)
#endif
{
	ProfileNodeManager::add(this);
}

void AutoProfileNode::release()
{
	if(m_bActive){
#ifdef UG_PROFILER_SHINY
		Shiny::ProfileManager::instance._endCurNode();
		PROFILE_LOG_CALL_END();
#endif
#ifdef UG_PROFILER_SCALASCA
		EPIK_USER_END(m_pName);
#endif
#ifdef UG_PROFILER_VAMPIR
		VT_USER_END((char*)m_pName);
#endif
#ifdef UG_PROFILER_SCOREP
		SCOREP_USER_REGION_END(m_pHandle);
#endif
		m_bActive = false;
	}
}

AutoProfileNode::~AutoProfileNode()
{
	if(m_bActive){
		ProfileNodeManager::release_latest();
	}
}

#endif
