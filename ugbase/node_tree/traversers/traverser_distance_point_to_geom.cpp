//	distance_point_to_geom.h
//	created by Sebastian Reiter y07 m12 d5
//	s.b.reiter@googlemail.com

#include <iostream>
#include "traverser_distance_point_to_geom.h"
#include "../node_tree.h"

using namespace std;

namespace ug{
namespace node_tree
{

Traverser_DistancePointToGeom::Traverser_DistancePointToGeom()
{
	m_closestEdgeID = -1;
}

Traverser_DistancePointToGeom::~Traverser_DistancePointToGeom()
{
}

number Traverser_DistancePointToGeom::
get_distance_point_to_geom(vector3& point, SPNode& nodeGraph)
{
	m_searchState = SEARCHING;
	m_distance = -1.0;
	m_closestEdgeID = -1;
	m_point = point;
	
	Traverser::apply(nodeGraph);

//	we found one
//	call the found-function
	if(m_searchState == GOT_ONE){
		closest_edge_found(m_closestEdgePointIndices[0],
							m_closestEdgePointIndices[1], 
							m_closestEdgeID, m_closestRootNode);
	}
	
	return m_distance;
}

int Traverser_DistancePointToGeom::get_closest_element_id()
{
	return m_closestEdgeID;
}

void Traverser_DistancePointToGeom::handle_group(GroupNode* group)
{	
	switch(m_searchState)
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
							Traverser::traverse_object(group->get_child(j).get_impl());
					}			
				}
			}
			break;

		case GOT_ONE:
			{
			//	traverse all children
				Traverser::handle_group(group);
			}
			break;
	}
}

void Traverser_DistancePointToGeom::handle_boxed_group(BoxedGroupNode* boxedGroup)
{
	switch(m_searchState)
	{
		case SEARCHING:
			//	check if the point lies inside the box. If so traverse the group.
			if(BoxBoundProbe(m_point, boxedGroup->min_corner(),
							boxedGroup->max_corner()))
			{
				Traverser::handle_boxed_group(boxedGroup);
			}
			break;

		case FORCE_FIND:
			Traverser::handle_boxed_group(boxedGroup);
			break;

		case GOT_ONE:
			if(BoxBoxIntersection(boxedGroup->min_corner(), boxedGroup->max_corner(),
									m_boxMin, m_boxMax))
			{
				Traverser::handle_boxed_group(boxedGroup);
			}
			break;
	}
}

void Traverser_DistancePointToGeom::
handle_collision_tree_root(CollisionTreeRootNode* collisionTreeRoot)
{
//	put the rootNode on top of the stack. All CollisionEdgeNodes children will index it's points
	m_stackRootNodes.push(collisionTreeRoot);

//	traverse the node
	Traverser::handle_group(collisionTreeRoot);
//	if we didn't find an edge force it
	if(m_searchState == SEARCHING)
	{
		m_searchState = FORCE_FIND;
		Traverser::handle_collision_tree_root(collisionTreeRoot);
	}

//	pop the rootNode from the stack
	m_stackRootNodes.pop();
}

void Traverser_DistancePointToGeom::
handle_collision_edges(CollisionEdgesNode* collisionEdges)
{
	CollisionTreeRootNode* root = m_stackRootNodes.top();
	const vector3* pPoints = root->get_points();
	int numIndices = collisionEdges->num_edges() * 2;
	const int* indices = collisionEdges->get_edges();

	for(int i = 0; i < numIndices; i+=2)
	{
		const vector3& p1 = pPoints[indices[i]];
		const vector3& p2 = pPoints[indices[i+1]];
		
		number distance = DistancePointToLine(m_point, p1, p2);

		if((m_searchState != GOT_ONE) || (distance < m_distance))
		{
			m_searchState = GOT_ONE;
			m_distance = distance;
			m_closestEdgePointIndices[0] = indices[i];
			m_closestEdgePointIndices[1] = indices[i+1];
			m_closestEdgeID = collisionEdges->get_edge_id(i/2);
			m_closestRootNode = root;
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

void Traverser_DistancePointToGeom::
closest_edge_found(int vrtInd1, int vrtInd2, int edgeID,
					CollisionTreeRootNode* collisionTree)
{
}

}//	end of namespace node_tree
}//	end of namespace ug
