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

#ifndef __H__UG__NODE_TREE__OBJECT__
#define __H__UG__NODE_TREE__OBJECT__

namespace ug {
namespace node_tree {


////////////////////////////////////////////////////////////////////////
//	Object-Codes
///	ids associated with objects. Values shouldn't be unnessesarily high.
enum ObjectCode
{
	OC_INVALID = 0xFFFFFFFF,
	OC_OBJECT = 0,
//	nodes
	OC_NODE = 1,
	OC_GROUP_NODE,
	OC_BOXED_GROUP_NODE,
	OC_COLLISION_TREE_ROOT_NODE,
	OC_COLLISION_EDGES_NODE,
	OC_COLLISION_TRIANGLES_NODE,

	OC_NODES_END,

//	custom objects
	OC_CUSTOM_OBJECT
};


////////////////////////////////////////////////////////////////////////
//	Object
///	An Object serves as the base-class for most of the polymorphic node-tree objects
/**
 * Each derivative should have its own ObjectCode, with which it can
 * be uniquely identified.
 */
class Object
{
	public:
		virtual ~Object()  = default;

		[[nodiscard]] inline unsigned int getObjectCode() const {return m_objectCode;}

	protected:
		Object() = default;
		Object(const Object& obj) {};

	protected:
		unsigned int m_objectCode;

};

}//	end of namespace node_tree
}//	end of namespace ug

#endif
