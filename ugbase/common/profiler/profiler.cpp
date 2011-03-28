// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 28.03.2011 (m,d,y)
 
#include "profiler.h"

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
		node->deactivate();
		Shiny::ProfileManager::instance._endCurNode();
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




AutoProfileNode::AutoProfileNode() : m_bActive(true)
{
	ProfileNodeManager::add(this);
}

AutoProfileNode::~AutoProfileNode()
{
	if(m_bActive){
		ProfileNodeManager::release_latest();
	}
}

