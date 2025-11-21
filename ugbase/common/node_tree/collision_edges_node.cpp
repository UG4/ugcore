/*
 * Copyright (c) 2007-2015:  Sebastian Reiter
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
void CollisionEdgesNode::add_edge(int ind1, int ind2)
{
	if(m_bEdgeIDsSupplied)
	{
		add_edge(ind1, ind2, CollisionElementID());
	}
	else
	{
		m_vEdges.push_back(ind1);
		m_vEdges.push_back(ind2);
	}
}

////////////////////////////////////////////////////////////////////////
void CollisionEdgesNode::
add_edge(int ind1, int ind2, CollisionElementID edgeID)
{
	if(!m_bEdgeIDsSupplied)
	{
		m_bEdgeIDsSupplied = true;
		if(num_edges() > 0)
			m_vEdgeIDs.resize(num_edges());
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
		m_vEdgeIDs.resize(maxEdge + numEdges);
	
	int startEdge = m_vEdges.size();
//	resize m_vEdges
	m_vEdges.resize(m_vEdges.size() + numEdges * 2);
	
	int upperBorder = numEdges*2;
	for(int i = 0; i < upperBorder; ++i)
		m_vEdges[startEdge + i] = pIndices[i];
}

////////////////////////////////////////////////////////////////////////
void CollisionEdgesNode::
add_edges(int* pIndices, CollisionElementID* pEdgeIDs, int numEdges)
{
	int maxEdge = m_vEdges.size() / 2;
	
	if(m_bEdgeIDsSupplied)
		m_vEdgeIDs.resize(maxEdge + numEdges);
	else
	{
		m_vEdgeIDs.resize(maxEdge + numEdges);
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
void CollisionEdgesNode::
set_edge_id(int edgeIndex, CollisionElementID edgeID)
{
	if(!m_bEdgeIDsSupplied)
	{
		m_bEdgeIDsSupplied = true;
		if(num_edges() > 0)
			m_vEdgeIDs.resize(num_edges());
	}

	m_vEdgeIDs[edgeIndex] = edgeID;
}

////////////////////////////////////////////////////////////////////////
CollisionElementID CollisionEdgesNode::
get_edge_id(int edgeIndex)
{
	if(m_bEdgeIDsSupplied)
		return m_vEdgeIDs[edgeIndex];

	return CollisionElementID();
}

}//	end of namespace node_tree
}//	end of namespace ug

