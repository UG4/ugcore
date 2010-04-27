//	distance_point_to_geom.h
//	created by Sebastian Reiter y07 m12 d5
//	s.b.reiter@googlemail.com

#include <iostream>
#include "traverser_project_point.h"
#include "../node_tree.h"

using namespace std;

namespace ug{
namespace node_tree
{

Traverser_ProjectPoint::Traverser_ProjectPoint()
{
	m_closestElemID = -1;
	m_closestElemType = 0;
}

Traverser_ProjectPoint::~Traverser_ProjectPoint()
{
}

bool Traverser_ProjectPoint::
project(const vector3& point, SPNode nodeGraph)
{
	m_searchState = SEARCHING;
	m_distance = -1.0;
	m_closestElemID = -1;
	m_closestElemType = 0;
	m_point = point;
	
	Traverser_CollisionTree::apply(nodeGraph);

	return m_searchState == GOT_ONE;
}

CollisionElementID Traverser_ProjectPoint::
get_closest_element_id()
{
	return m_closestElemID;
}

void Traverser_ProjectPoint::handle_group(GroupNode* group)
{	
	SearchState tSearchState = m_searchState;
	switch(tSearchState)
	{
		case SEARCHING:
		case FORCE_FIND:
			{
			//	traverse each subnode
				int numChildren = group->num_children();
				int i;

				for(i = 0; i < numChildren; ++i)
				{
					traverse_object(group->get_child(i).get_impl());
				//	check if the gotOne state changed during the traversal of the sub-node
					if(m_searchState == GOT_ONE){
					//	it did. break the loop and revisit all children of the group.
						break;
					}
				}

				if(m_searchState == GOT_ONE)
				{
				//	the got one state was just enabled.
				//	We have to recheck all children with m_gotOne enabled.
				//	ignore child i (since this has been checked already)
					for(int j = 0; j < numChildren; ++j)
					{
						if(j != i)
							Traverser_CollisionTree::traverse_object(group->get_child(j).get_impl());
					}			
				}
			}
			break;

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
			//	check if the point lies inside the box. If so traverse the group.
			if(BoxBoundProbe(m_point, boxedGroup->min_corner(),
							boxedGroup->max_corner()))
			{
				Traverser_CollisionTree::handle_boxed_group(boxedGroup);
			}
			break;

		case FORCE_FIND:
			Traverser_CollisionTree::handle_boxed_group(boxedGroup);
			break;

		case GOT_ONE:
			if(BoxBoxIntersection(boxedGroup->min_corner(), boxedGroup->max_corner(),
									m_boxMin, m_boxMax))
			{
				Traverser_CollisionTree::handle_boxed_group(boxedGroup);
			}
			break;
	}
}

void Traverser_ProjectPoint::
handle_collision_tree_root(CollisionTreeRootNode* colTreeRootNode)
{
//	put the rootNode on top of the stack.
//	All CollisionEdgeNodes children will index it's points
	m_stackRootNodes.push(colTreeRootNode);

//	traverse the node
	Traverser_CollisionTree::handle_group(colTreeRootNode);
//	if we didn't find an edge force it
	if(m_searchState == SEARCHING)
	{
		m_searchState = FORCE_FIND;
		Traverser_CollisionTree::handle_collision_tree_root(colTreeRootNode);
	}

//	pop the rootNode from the stack
	m_stackRootNodes.pop();
}

void Traverser_ProjectPoint::
handle_collision_edges(CollisionEdgesNode* colEdgesNode)
{
	CollisionTreeRootNode* root = m_stackRootNodes.top();
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
			m_boxMin.x = m_point.x - distance;			
			m_boxMin.y = m_point.y - distance;
			m_boxMin.z = m_point.z - distance;
			m_boxMax.x = m_point.x + distance;
			m_boxMax.y = m_point.y + distance;
			m_boxMax.z = m_point.z + distance;
		}
	}
}

void Traverser_ProjectPoint::
handle_collision_triangles(CollisionTrianglesNode* colTrisNode)
{
	CollisionTreeRootNode* root = m_stackRootNodes.top();
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
			m_boxMin.x = m_point.x - distance;			
			m_boxMin.y = m_point.y - distance;
			m_boxMin.z = m_point.z - distance;
			m_boxMax.x = m_point.x + distance;
			m_boxMax.y = m_point.y + distance;
			m_boxMax.z = m_point.z + distance;
		}
	}
}

}//	end of namespace node_tree
}//	end of namespace ug
