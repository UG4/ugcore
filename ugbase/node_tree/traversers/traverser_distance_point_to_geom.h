//	distance_point_to_geom.h
//	created by Sebastian Reiter y07 m12 d5
//	s.b.reiter@googlemail.com

#ifndef __H__UG__NODE_TREE__TRAVERSER_DISTANCE_POINT_TO_GEOM__
#define __H__UG__NODE_TREE__TRAVERSER_DISTANCE_POINT_TO_GEOM__

#include <stack>
#include "../traverser.h"
#include "common/math/ugmath.h"

namespace ug{
namespace node_tree
{
////////////////////////////////////////////////////////////////////////
///	traverses a node-tree and projects a given point to the contained geometry.
class Traverser_DistancePointToGeom : protected Traverser
{
	public:
		Traverser_DistancePointToGeom();
		virtual ~Traverser_DistancePointToGeom();

		virtual number get_distance_point_to_geom(vector3& point,
												  SPNode& nodeGraph);

	/** after the distance of a point to the geometry has been determined,
	 *	this funtion returns the id of the element closest to the point.*/
		virtual int get_closest_element_id();
	//	as soon as triangles are supported a function has to be added
	//	that returns the type of the closest element!
		
	protected:
		virtual void handle_group(GroupNode* group);
		virtual void handle_boxed_group(BoxedGroupNode* boxedGroup);
		virtual void handle_collision_tree_root(CollisionTreeRootNode* collisionTreeRoot);
		virtual void handle_collision_edges(CollisionEdgesNode* collisionEdges);

	/**	this funtion is called for each completly processed CollisionTreeRootNode
	 *	(if the node contained edges).*/
		virtual void closest_edge_found(int vrtInd1, int vrtInd2, int edgeID,
										CollisionTreeRootNode* collisionTree);

	protected:
		enum SearchState
		{
			SEARCHING = 0,
			FORCE_FIND,
			GOT_ONE
		};

		vector3	m_point;
		number	m_distance;

	private:
		SearchState			m_searchState;

		vector3	m_boxMin;
		vector3	m_boxMax;

		int					m_closestEdgePointIndices[2];
		int					m_closestEdgeID;
		CollisionTreeRootNode*	m_closestRootNode;

		std::stack<CollisionTreeRootNode*>	m_stackRootNodes;
};

}//	end of namespace node_tree
}//	end of namespace ug

#endif
