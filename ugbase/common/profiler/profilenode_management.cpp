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
ProfileNodeManager()	{}

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
