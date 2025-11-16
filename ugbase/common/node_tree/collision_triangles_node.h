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

#ifndef __H__UG__NODE_TREE__COLLSISION_TRIANGLES_NODE__
#define __H__UG__NODE_TREE__COLLSISION_TRIANGLES_NODE__

#include <vector>
#include "node.h"
#include "collision_element_info.h"

namespace ug{
namespace node_tree
{
class CollisionTrianglesNode;

////////////////////////////////////////////////////////////////////////
///	the smartpointer used to encapsulate the node
using SPCollisionTrianglesNode = SmartPtr<CollisionTrianglesNode>;

////////////////////////////////////////////////////////////////////////
//	CollisionTrianglesNode
///	holds index tuples defining triangles.
/**
 * The index-tuples refer to the next CollisionTreeRootNode, which is
 * higher in the hierarchy.
 *
 * An identifier can be stored with each triangle - either an
 * integer-value or a void-pointer.
 * Normally appears as a subordinate of CollisionTreeRootNode.
 */
class CollisionTrianglesNode : public Node
{
	public:
		static SPCollisionTrianglesNode create();

		virtual ~CollisionTrianglesNode();

		virtual void add_triangle(int ind1, int ind2, int ind3);
		virtual void add_triangle(int ind1, int ind2, int ind3,
								  CollisionElementID triID);
		
	/// pIndices has to be of size numTris*3
		virtual void add_triangles(int* pIndices, size_t numTris);
		
	///	pIndices and pTriIDs have to be of size numTris*3
		virtual void add_triangles(int* pIndices,
								   CollisionElementID* pTriIDs,
								   size_t numTris);
		
		virtual size_t num_triangles() const;
		
		virtual void get_triangle(size_t index, int& ind1Out,
								  int& ind2Out, int& ind3Out) const;

		virtual const int* get_triangles() const;

		virtual void set_triangle_id(size_t triInd,
									 CollisionElementID triID);
		
	/// if no identifier has been set for an edge -1 is returned.
		virtual CollisionElementID get_triangle_id(size_t triInd);

	protected:
		CollisionTrianglesNode();

	protected:
		using IndexVec = std::vector<int>;
		using IDVec = std::vector<CollisionElementID>;

		IndexVec	m_vTris;
		IDVec		m_vTriIDs;
		
		bool m_bTriangleIDsSupplied;
};

}//	end of namespace node_tree
}//	end of namespace ug

#endif
