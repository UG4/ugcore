//	collision_edges_node.cpp
//	created by Sebastian Reiter y07 m12 d5
//	s.b.reiter@googlemail.com

#include "collision_edges_node.h"

namespace ug{
namespace node_tree
{

////////////////////////////////////////////////////////////////////////
SPCollisionEdgesNode CollisionEdgesNode::create()
{
	CollisionEdgesNode* node = new CollisionEdgesNode;
	node->m_objectCode = OC_COLLISION_EDGES_NODE;
	return SPCollisionEdgesNode(node);
}

////////////////////////////////////////////////////////////////////////
CollisionEdgesNode::CollisionEdgesNode()
{
	m_bEdgeIDsSupplied = false;
}

////////////////////////////////////////////////////////////////////////
CollisionEdgesNode::~CollisionEdgesNode()
{
}

////////////////////////////////////////////////////////////////////////
void CollisionEdgesNode::add_edge(int ind1, int ind2)
{
	if(m_bEdgeIDsSupplied)
	{
		add_edge(ind1, ind2, -1);
	}
	else
	{
		m_vEdges.push_back(ind1);
		m_vEdges.push_back(ind2);
	}
}

////////////////////////////////////////////////////////////////////////
void CollisionEdgesNode::add_edge(int ind1, int ind2, int edgeID)
{
	if(!m_bEdgeIDsSupplied)
	{
		m_bEdgeIDsSupplied = true;
		if(num_edges() > 0)
			m_vEdgeIDs.resize(num_edges(), -1);
	}

	m_vEdges.push_back(ind1);
	m_vEdges.push_back(ind2);
	m_vEdgeIDs.push_back(edgeID);
}

////////////////////////////////////////////////////////////////////////
void CollisionEdgesNode::add_edges(int* pIndices, int numEdges)
{
	int maxEdge = m_vEdges.size() / 2;
	
	if(m_bEdgeIDsSupplied)
		m_vEdgeIDs.resize(maxEdge + numEdges, -1);
	
	int startEdge = m_vEdges.size();
//	resize m_vEdges
	m_vEdges.resize(m_vEdges.size() + numEdges * 2);
	
	int upperBorder = numEdges*2;
	for(int i = 0; i < upperBorder; ++i)
		m_vEdges[startEdge + i] = pIndices[i];
}

////////////////////////////////////////////////////////////////////////
void CollisionEdgesNode::add_edges(int* pIndices, int* pEdgeIDs, int numEdges)
{
	int maxEdge = m_vEdges.size() / 2;
	
	if(m_bEdgeIDsSupplied)
		m_vEdgeIDs.resize(maxEdge + numEdges);
	else
	{
		m_vEdgeIDs.resize(maxEdge + numEdges, -1);
		m_bEdgeIDsSupplied = true;
	}

	
//	fill the ids
	{
		for(int i = 0; i <numEdges; ++i)
			m_vEdgeIDs[maxEdge + i] = pEdgeIDs[i];
	}	
	
	int startEdge = m_vEdges.size();
//	resize m_vEdges
	m_vEdges.resize(m_vEdges.size() + numEdges * 2);
	
//	fill m_vEdges
	{
		int upperBorder = numEdges*2;
		for(int i = 0; i < upperBorder; ++i)
			m_vEdges[startEdge + i] = pIndices[i];
	}

}

////////////////////////////////////////////////////////////////////////
int CollisionEdgesNode::num_edges() const
{
	return(m_vEdges.size() / 2);
}

////////////////////////////////////////////////////////////////////////
void CollisionEdgesNode::get_edge(int index, int& ind1Out, int& ind2Out) const
{
	index *= 2;
	ind1Out = m_vEdges[index++];
	ind2Out = m_vEdges[index];
}

////////////////////////////////////////////////////////////////////////
const int* CollisionEdgesNode::get_edges() const
{
	return &m_vEdges.front();
}

////////////////////////////////////////////////////////////////////////
void CollisionEdgesNode::set_edge_id(int edgeIndex, int edgeID)
{
	if(!m_bEdgeIDsSupplied)
	{
		m_bEdgeIDsSupplied = true;
		if(num_edges() > 0)
			m_vEdgeIDs.resize(num_edges(), -1);
	}

	m_vEdgeIDs[edgeIndex] = edgeID;
}

////////////////////////////////////////////////////////////////////////
int CollisionEdgesNode::get_edge_id(int edgeIndex)
{
	if(m_bEdgeIDsSupplied)
		return m_vEdgeIDs[edgeIndex];

	return -1;
}

}//	end of namespace node_tree
}//	end of namespace ug

