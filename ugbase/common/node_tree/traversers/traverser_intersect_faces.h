/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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
