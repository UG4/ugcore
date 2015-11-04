#ifndef __H__UG__NODE_TREE__INTERSECT_FACES__
#define __H__UG__NODE_TREE__INTERSECT_FACES__

#include <vector>
#include "common/math/ugmath.h"
#include "traverser_collision_tree.h"

namespace ug{
namespace node_tree
{
////////////////////////////////////////////////////////////////////////////////
///	traverses a node-tree and intersect a given face with the contained geometry.
/**	Currently only triangles are supported.
 *
 * \todo	A list with local coordinates of intersections should be created
 * 			during intersect_tri. This list should then be available to the user.
 */
class Traverser_IntersectFaces : protected Traverser_CollisionTree
{
	public:
		Traverser_IntersectFaces();
		virtual ~Traverser_IntersectFaces();

	///	intersects the given triangle with all faces in the given nodeGraph.
	/**	returns true if an intersection was found, false if not.
	 * After each run the intersected faces can be accessed using
	 * get_intersected_element_ids().*/
		virtual bool intersect_tri(const vector3& v0, const vector3& v1,
									const vector3& v2, SPNode nodeGraph);

	//todo:	Add intersect_quad(...)

	///	adds an element to the ignore list
	/**	Make sure that the ignore list won't get too big, since for each
	 * it has to be traversed for each triangle that intersects.
	 * Use clear_ignore_list to clear the list.*/
		void ignore_element(const CollisionElementID& elemID);

	///	clears the ignore list
		void clear_ignore_list();

	/** after the intersection with the geometry has been performed,
	 *	this function returns the ids of the intersected elements.*/
		const std::vector<CollisionElementID>& get_intersected_element_ids() const;
		
	protected:
		virtual void handle_boxed_group(BoxedGroupNode* boxedGroup);
		virtual void handle_collision_triangles(CollisionTrianglesNode* colTrisNode);
		
	private:
	//	the element which shall be checked
		vector3	m_vrts[4];
		int		m_numVrts;

	//	the intersecting elements will be stored here
		std::vector<CollisionElementID>	m_intersectedElementIDs;

	//	an intersection is only recorded if the intersecting element is not
	//	contained in the ignore list.
		std::vector<CollisionElementID>	m_ignoreList;
};

}//	end of namespace node_tree
}//	end of namespace ug

#endif
