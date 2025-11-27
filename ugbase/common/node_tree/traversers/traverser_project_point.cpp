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

#include <iostream>
#include "traverser_project_point.h"
#include "../node_tree.h"
#include "common/log.h"
#include "common/profiler/profiler.h"

using namespace std;

namespace ug{
namespace node_tree
{

/*
force_find muss überarbeitet werden:
Problem: kleinste box in der der knoten liegt hält keine dreiecke.
Ansatz: Falls nach dem Abarbeiten der nächst kleineren Box der Status
	noch auf SEARCHING steht, wird er auf FORCE_FIND gesetzt.
	War der Status auf SEARCHING, wurde aber von Kindknoten oder dem
	Knoten selbst auf FORCE_FIND gesetzt, so wird über alle Nachbarn
	der Box iteriert und mit ihnen ein Test durchgeführt.
	Ist der Status zu Begin der Evaluierung eines Knotens auf FORCE_FIND,
	so werden gleich alle Kind-Knoten evaluiert.
*/

Traverser_ProjectPoint::Traverser_ProjectPoint()
{
	m_closestElemID = -1;
	m_closestElemType = 0;
}


bool Traverser_ProjectPoint::
project(const vector3& point, SPNode nodeGraph,
		vector3* pPointNormal)
{
	m_searchState = SEARCHING;
	m_distance = -1.0;
	m_closestElemID = -1;
	m_closestElemType = 0;
	m_point = point;
	
	if(pPointNormal){
		m_checkNormals = true;
		m_pointNormal = *pPointNormal;
	}
	else
		m_checkNormals= false;
	
	apply(nodeGraph);

	return m_searchState == GOT_ONE;
}

CollisionElementID Traverser_ProjectPoint::
get_closest_element_id() const {
	return m_closestElemID;
}

void Traverser_ProjectPoint::handle_group(GroupNode* group)
{
	SearchState tSearchState = m_searchState;
	switch(tSearchState)
	{
		case SEARCHING:
			{
			//	traverse each subnode
				int numChildren = group->num_children();

				for(int i = 0; i < numChildren; ++i)
				{
					traverse_object(group->get_child(i).get());
				
				//	check whether the state changed during traversal of
				//	the i-th subnode
					if(m_searchState != SEARCHING){
					//	the state changed. we have to traverse the neighbours
					//	of this node again, since other boxed-groups may be
					//	entered now.
					//	ignore child i (since this has been checked already)
						for(int j = 0; j < numChildren; ++j)
						{
							if(j != i)
								traverse_object(group->get_child(j).get());
						}
					//	all children have been traversed. We're done in here.
						break;
					}
				}
			
			//	At this point we have traversed all child nodes of
			//	the node. If we're still in SEARCHING mode at this
			//	point, we have to force a find and re-traverse the
			//	children.
				if(m_searchState == SEARCHING){
					m_searchState = FORCE_FIND;
					Traverser_CollisionTree::handle_group(group);
				}
			}
			break;
			
		case FORCE_FIND:
		case GOT_ONE:
			{
			//	traverse all children
				Traverser_CollisionTree::handle_group(group);
			}
			break;
	}
}

void Traverser_ProjectPoint::handle_boxed_group(BoxedGroupNode* boxedGroup)
{
	SearchState tSearchState = m_searchState;
	switch(tSearchState)
	{
		case SEARCHING:
		//	check whether the point lies inside the box.
		//	If so traverse the group (all children).
			if(BoxBoundProbe(m_point, boxedGroup->min_corner(),
							 boxedGroup->max_corner()))
			{
				handle_group(boxedGroup);
			}
			break;

		case FORCE_FIND:
			handle_group(boxedGroup);
			break;

		case GOT_ONE:
		//	check whether the box around the point being projected, which
		//	contains all points which could possibly be closer than the
		//	temporarily closest point, intersects the boxedGroups bounding box.
			if(BoxBoxIntersection(boxedGroup->min_corner(),
								  boxedGroup->max_corner(),
								  m_boxMin, m_boxMax))
			{
				handle_group(boxedGroup);
			}
			break;
	}
}

void Traverser_ProjectPoint::
handle_collision_edges(CollisionEdgesNode* colEdgesNode)
{
	CollisionTreeRootNode* root = get_current_root_node();
	const vector3* pPoints = root->get_points();
	int numIndices = colEdgesNode->num_edges() * 2;
	const int* indices = colEdgesNode->get_edges();

	for(int i = 0; i < numIndices; i+=2)
	{
		const vector3& p1 = pPoints[indices[i]];
		const vector3& p2 = pPoints[indices[i+1]];
		
		number t;
		number distance = DistancePointToLine(t, m_point, p1, p2);
		
		if((m_searchState != GOT_ONE) || (distance < m_distance))
		{
			m_searchState = GOT_ONE;
			m_distance = distance;
			m_closestElemIndices[0] = indices[i];
			m_closestElemIndices[1] = indices[i+1];
			m_closestElemID = colEdgesNode->get_edge_id(i/2);
			m_closestElemType = 2;
			m_closestRootNode = root;
			VecScaleAdd(m_closestPoint, 1. - t, p1, t, p2);
			
		//	reset the box
			m_boxMin.x() = m_point.x() - distance * 1.01;			
			m_boxMin.y() = m_point.y() - distance * 1.01;
			m_boxMin.z() = m_point.z() - distance * 1.01;
			m_boxMax.x() = m_point.x() + distance * 1.01;
			m_boxMax.y() = m_point.y() + distance * 1.01;
			m_boxMax.z() = m_point.z() + distance * 1.01;
		}
	}
}

void Traverser_ProjectPoint::
handle_collision_triangles(CollisionTrianglesNode* colTrisNode)
{
	CollisionTreeRootNode* root = get_current_root_node();
	const vector3* pPoints = root->get_points();
	int numIndices = colTrisNode->num_triangles() * 3;
	const int* indices = colTrisNode->get_triangles();

	for(int i = 0; i < numIndices; i+=3)
	{
		const vector3& p1 = pPoints[indices[i]];
		const vector3& p2 = pPoints[indices[i+1]];
		const vector3& p3 = pPoints[indices[i+2]];
		
		vector3 n;
		CalculateTriangleNormalNoNormalize(n, p1, p2, p3);
		
	//	if normal-check is enabled, we have to make sure, that the points
	//	normal points into the same direction as the triangles normal.
		if(m_checkNormals){
			if(VecDot(n, m_pointNormal) <= 0)
				continue;
		}
		
		number bc1, bc2;
		vector3 vTmp;
		number distance = DistancePointToTriangle(vTmp, bc1, bc2,
													m_point, p1, p2, p3, n);
				
		if((m_searchState != GOT_ONE) || (distance < m_distance))
		{
			m_searchState = GOT_ONE;
			m_distance = distance;
			m_closestElemIndices[0] = indices[i];
			m_closestElemIndices[1] = indices[i+1];
			m_closestElemIndices[2] = indices[i+2];
			m_closestElemID = colTrisNode->get_triangle_id(i/3);
			m_closestElemType = 3;
			m_closestRootNode = root;
			m_closestPoint = vTmp;
			
		//	reset the box
			m_boxMin.x() = m_point.x() - distance * 1.01;			
			m_boxMin.y() = m_point.y() - distance * 1.01;
			m_boxMin.z() = m_point.z() - distance * 1.01;
			m_boxMax.x() = m_point.x() + distance * 1.01;
			m_boxMax.y() = m_point.y() + distance * 1.01;
			m_boxMax.z() = m_point.z() + distance * 1.01;
		}
	}
}

}//	end of namespace node_tree
}//	end of namespace ug
