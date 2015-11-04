#ifndef __H__UG__NODE_TREE__TRAVERSER_PROJECT_POINT__
#define __H__UG__NODE_TREE__TRAVERSER_PROJECT_POINT__

#include "common/math/ugmath.h"
#include "traverser_collision_tree.h"

namespace ug{
namespace node_tree
{
////////////////////////////////////////////////////////////////////////
///	traverses a node-tree and projects a given point to the contained geometry.
class Traverser_ProjectPoint : protected Traverser_CollisionTree
{
	public:
		Traverser_ProjectPoint();
		virtual ~Traverser_ProjectPoint();

		virtual bool project(const vector3& point, SPNode nodeGraph,
							vector3* pPointNormal = NULL);

	/** after the distance of a point to the geometry has been determined,
	 *	this funtion returns the id of the element closest to the point.*/
		CollisionElementID get_closest_element_id();
		
	///	returns the type of the element to which the closest distance was found.
	/**	0: invalid,
	 *	2: edge,
	 *	3: triangle*/
	 	int get_closest_element_type();
		
		inline number get_distance()			{return m_distance;}
		
		inline vector3 get_closest_point()		{return m_closestPoint;}
		
	protected:
		virtual void handle_group(GroupNode* group);
		virtual void handle_boxed_group(BoxedGroupNode* boxedGroup);
		virtual void handle_collision_edges(CollisionEdgesNode* colTrisNode);
		virtual void handle_collision_triangles(CollisionTrianglesNode* colTrisNode);
		
	protected:
		enum SearchState
		{
			SEARCHING = 0,
			FORCE_FIND,
			GOT_ONE
		};

		vector3	m_point;
		vector3 m_closestPoint;
		number	m_distance;

	private:
		SearchState			m_searchState;
		bool 	m_checkNormals;
		vector3	m_pointNormal;
		
		vector3	m_boxMin;
		vector3	m_boxMax;

		int						m_closestElemIndices[3];
		int						m_closestElemType;
		CollisionElementID		m_closestElemID;
		CollisionTreeRootNode*	m_closestRootNode;
};

}//	end of namespace node_tree
}//	end of namespace ug

#endif
