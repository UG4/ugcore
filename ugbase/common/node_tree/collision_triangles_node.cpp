/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#include "collision_triangles_node.h"

namespace ug {
namespace node_tree {

////////////////////////////////////////////////////////////////////////
SPCollisionTrianglesNode CollisionTrianglesNode::create()
{
	auto node = new CollisionTrianglesNode;
	node->m_objectCode = OC_COLLISION_TRIANGLES_NODE;
	return SPCollisionTrianglesNode(node);
}

////////////////////////////////////////////////////////////////////////
CollisionTrianglesNode::CollisionTrianglesNode()
{
	m_bTriangleIDsSupplied = false;
}


////////////////////////////////////////////////////////////////////////
void CollisionTrianglesNode::add_triangle(int ind1, int ind2, int ind3)
{
	if(m_bTriangleIDsSupplied)
	{
		add_triangle(ind1, ind2, ind3);
	}
	else
	{
		m_vTris.push_back(ind1);
		m_vTris.push_back(ind2);
		m_vTris.push_back(ind3);
	}
}

////////////////////////////////////////////////////////////////////////
void CollisionTrianglesNode::add_triangle(int ind1, int ind2, int ind3,
										  CollisionElementID triID)
{
	m_vTris.push_back(ind1);
	m_vTris.push_back(ind2);
	m_vTris.push_back(ind3);
	
	if(!m_bTriangleIDsSupplied)
	{
		m_bTriangleIDsSupplied = true;
		if(num_triangles() > 0)
			m_vTriIDs.resize(num_triangles());
	}
	
	m_vTriIDs.push_back(triID);
}

////////////////////////////////////////////////////////////////////////
void CollisionTrianglesNode::add_triangles(int* pIndices, size_t numTris)
{
	size_t maxTri = m_vTris.size() / 3;
	
	if(m_bTriangleIDsSupplied)
		m_vTriIDs.resize(maxTri + numTris);
	
	size_t startInd = m_vTris.size();
	
//	resize m_vTris
	m_vTris.resize(m_vTris.size() + numTris * 3);
	
	size_t numInds = numTris * 3;
	for(size_t i = 0; i < numInds; ++i)
		m_vTris[startInd + i] = pIndices[i];
}

////////////////////////////////////////////////////////////////////////
void CollisionTrianglesNode::
add_triangles(int* pIndices, CollisionElementID* pTriIDs, size_t numTris)
{
	size_t maxTri = m_vTris.size() / 3;
	
	if(m_bTriangleIDsSupplied)
		m_vTriIDs.resize(maxTri + numTris);
	else
	{
		m_vTriIDs.resize(maxTri + numTris);
		m_bTriangleIDsSupplied = true;
	}

	
//	fill the ids
	{
		for(size_t i = 0; i < numTris; ++i)
			m_vTriIDs[maxTri + i] = pTriIDs[i];
	}	
	
	size_t startInd = m_vTris.size();
	
//	resize m_vTris
	m_vTris.resize(m_vTris.size() + numTris * 3);
	
//	fill m_vEdges
	{
		size_t numInds = numTris * 3;
		for(size_t i = 0; i < numInds; ++i)
			m_vTris[startInd + i] = pIndices[i];
	}

}

////////////////////////////////////////////////////////////////////////
size_t CollisionTrianglesNode::num_triangles() const
{
	return(m_vTris.size() / 3);
}

////////////////////////////////////////////////////////////////////////
void CollisionTrianglesNode::
get_triangle(size_t index, int& ind1Out, int& ind2Out, int& ind3Out) const
{
	index *= 3;
	ind1Out = m_vTris[index++];
	ind2Out = m_vTris[index++];
	ind3Out = m_vTris[index];
}

////////////////////////////////////////////////////////////////////////
const int* CollisionTrianglesNode::get_triangles() const
{
	return &m_vTris.front();
}

////////////////////////////////////////////////////////////////////////
void CollisionTrianglesNode::
set_triangle_id(size_t triInd, CollisionElementID triID)
{
	if(!m_bTriangleIDsSupplied)
	{
		m_bTriangleIDsSupplied = true;
		if(num_triangles() > 0)
			m_vTriIDs.resize(num_triangles());
	}

	if(triInd >= m_vTriIDs.size())
		m_vTriIDs.resize(triInd + 1);
		
	m_vTriIDs[triInd] = triID;
}

////////////////////////////////////////////////////////////////////////
CollisionElementID
CollisionTrianglesNode::get_triangle_id(size_t triInd)
{
	if(m_bTriangleIDsSupplied)
		return m_vTriIDs[triInd];

	return CollisionElementID();
}

}//	end of namespace node_tree
}//	end of namespace ug

