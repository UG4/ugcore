// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 28.03.2011 (m,d,y)
 
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




AutoProfileNode::AutoProfileNode(const char* name)
	: m_bActive(true), m_pName(name)
{
	ProfileNodeManager::add(this);
}

void AutoProfileNode::release()
{
	if(m_bActive){
#ifdef UG_PROFILER_SHINY
		Shiny::ProfileManager::instance._endCurNode();
#endif
#ifdef UG_PROFILER_SCALASCA
		EPIK_USER_END(m_pName);
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
