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

#ifndef __H__UG__NODE_TREE__TRAVERSER__
#define __H__UG__NODE_TREE__TRAVERSER__

#include <vector>
#include "object.h"
#include "node.h"

namespace ug{
namespace node_tree
{

class Traverser;
class GroupNode;
class BoxedGroupNode;

////////////////////////////////////////////////////////////////////////
//	Traverser
///	Derivates of a Traverser can be used to traverse a scenegraph.
/**
*/
class Traverser
{
	public:
		Traverser();
		virtual ~Traverser() = default;

	protected:
		void apply(SPNode& node);

		template<typename HandlerType>
		void register_handler_function(unsigned int oc, HandlerType func);

		void traverse_object(Object* obj);

		virtual void handle_group(GroupNode* group);
		virtual void handle_boxed_group(BoxedGroupNode* boxedGroup);

	private:
		[[nodiscard]] bool handler_function_registered(unsigned int oc) const;

	private:
		using HandlerFunc = void(Traverser::*)(Object* obj);
		std::vector<HandlerFunc>	m_vHandlerFuncs;
};


template<typename HandlerType>
void Traverser::register_handler_function(unsigned int oc, HandlerType func)
{
//	make sure that there is enough space
	if(oc >= m_vHandlerFuncs.size())
		m_vHandlerFuncs.resize(oc+1, nullptr);

	m_vHandlerFuncs[oc] = (HandlerFunc)func;
}


}//	end of namespace node_tree
}//	end of namespace ug

#endif
