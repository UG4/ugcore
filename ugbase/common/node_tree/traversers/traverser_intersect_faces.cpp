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

//#include <iostream>
#include "traverser_intersect_faces.h"
#include "../node_tree.h"
#include "common/log.h"
#include "common/profiler/profiler.h"

using namespace std;

namespace ug{
namespace node_tree
{


bool Traverser_IntersectFaces::
intersect_tri(const vector3& v0, const vector3& v1,
			  const vector3& v2, SPNode nodeGraph)
{
	m_vrts[0] = v0; m_vrts[1] = v1; m_vrts[2] = v2;
	m_numVrts = 3;
	
	m_intersectedElementIDs.clear();
	
	apply(nodeGraph);

	return !m_intersectedElementIDs.empty();
}

const std::vector<CollisionElementID>& Traverser_IntersectFaces::
get_intersected_element_ids() const
{
	return m_intersectedElementIDs;
}

void Traverser_IntersectFaces::
handle_boxed_group(BoxedGroupNode* boxedGroup)
{
//	check whether our element intersects the given box
//	todo: consider quadrilaterals
	if(TriangleBoxIntersection(m_vrts[0], m_vrts[1], m_vrts[2],
								boxedGroup->min_corner(),
								boxedGroup->max_corner()))
	{
		handle_group(boxedGroup);
	}
}

void Traverser_IntersectFaces::
handle_collision_triangles(CollisionTrianglesNode* colTrisNode)
{
	CollisionTreeRootNode* root = get_current_root_node();
	const vector3* pPoints = root->get_points();
	int numIndices = colTrisNode->num_triangles() * 3;
	const int* indices = colTrisNode->get_triangles();

//	iterate over all triangles of this node
	for(int i = 0; i < numIndices; i+=3)
	{
	//	todo: instead of checking the ignore list after intersection, it should
	//			probably be checked before. This depends on the size of the
	//			ignore list. Probably an ignore hash would be better.

	//	perform intersection
	//	todo: store local coordinates
		if(TriangleTriangleIntersection(m_vrts[0], m_vrts[1], m_vrts[2],
										pPoints[indices[i]],
										pPoints[indices[i+1]],
										pPoints[indices[i+2]]))
		{
		//	check whether the element is in the ignore list
			const CollisionElementID& id = colTrisNode->get_triangle_id(i/3);
			if(find(m_ignoreList.begin(), m_ignoreList.end(), id) == m_ignoreList.end()){
				m_intersectedElementIDs.push_back(id);
			}
		}
	}
}

void Traverser_IntersectFaces::
ignore_element(const CollisionElementID& elemID)
{
	m_ignoreList.push_back(elemID);
}

void Traverser_IntersectFaces::
clear_ignore_list()
{
	m_ignoreList.clear();
}

}//	end of namespace node_tree
}//	end of namespace ug
