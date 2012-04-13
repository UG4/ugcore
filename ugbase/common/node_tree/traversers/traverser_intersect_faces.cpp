//	distance_point_to_geom.h
//	created by Sebastian Reiter y07 m12 d5
//	s.b.reiter@googlemail.com

//#include <iostream>
#include "traverser_intersect_faces.h"
#include "../node_tree.h"
#include "common/log.h"
#include "common/profiler/profiler.h"

using namespace std;

namespace ug{
namespace node_tree
{

Traverser_IntersectFaces::Traverser_IntersectFaces()
{
}

Traverser_IntersectFaces::~Traverser_IntersectFaces()
{
}

bool Traverser_IntersectFaces::
intersect_tri(const vector3& v0, const vector3& v1,
			  const vector3& v2, SPNode nodeGraph)
{
	m_vrts[0] = v0; m_vrts[1] = v1; m_vrts[2] = v2;
	m_numVrts = 3;
	
	m_intersectedElementIDs.clear();
	
	Traverser_CollisionTree::apply(nodeGraph);

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
