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

#ifndef __H__UG__NODE_TREE__COLLSISION_EDGES_NODE__
#define __H__UG__NODE_TREE__COLLSISION_EDGES_NODE__

#include <vector>

#include "node.h"
#include "collision_element_info.h"

namespace ug {
namespace node_tree {

class CollisionEdgesNode;

////////////////////////////////////////////////////////////////////////
///	the smartpointer used to encapsulate the node
using SPCollisionEdgesNode = SmartPtr<CollisionEdgesNode>;

////////////////////////////////////////////////////////////////////////
//	CollisionEdgesNode
///	holds index pairs defining edges.
/**
 * An identifier can be stored with each edge.
 * Normally appears as a subordinate of CollisionTreeRootNode.
 */
class CollisionEdgesNode : public Node
{
	public:
		static SPCollisionEdgesNode create();

		~CollisionEdgesNode() override = default;

		virtual void add_edge(int ind1, int ind2);
		virtual void add_edge(int ind1, int ind2,
							  CollisionElementID edgeID);
		
	/// pIndices has to be of size numEdges*2
		virtual void add_edges(int* pIndices, int numEdges);
		
	///	pIndices has to be of size numEdges*2
		virtual void add_edges(int* pIndices,
							   CollisionElementID* pEdgeIDs,
							   int numEdges);
							   
		virtual int num_edges() const;
		
		virtual void get_edge(int index, int& ind1Out, int& ind2Out) const;
		
		virtual const int* get_edges() const;

		virtual void set_edge_id(int edgeIndex,
								 CollisionElementID edgeID);
		
	/// if no identifier has been set for an edge -1 is returned.
		virtual CollisionElementID get_edge_id(int edgeIndex);

	protected:
		CollisionEdgesNode();

	protected:
		using IndexVec = std::vector<int>;
		using IDVec = std::vector<CollisionElementID>;

		IndexVec	m_vEdges;
		IDVec		m_vEdgeIDs;

		bool m_bEdgeIDsSupplied;
};

}//	end of namespace node_tree
}//	end of namespace ug

#endif
