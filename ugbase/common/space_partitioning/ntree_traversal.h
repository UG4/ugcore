/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__ntree_traversal__
#define __H__UG__ntree_traversal__

namespace ug {

enum TraversalStates{
	DONT_TRAVERSE_CHILDREN,
	TRAVERSE_CHILDREN,
	ABORT_TRAVERSAL
};

template <typename tree_t, typename traverser_t>
void TraverseBreadthFirst(const tree_t& tree, traverser_t& traverser)
{
	std::vector<size_t> nodes;
	nodes.reserve(tree.num_nodes());
	nodes.push_back(0);
	int curInd = 0;

	traverser.begin_traversal(tree);

//	traverse upwards
	while(curInd < (int)nodes.size()){
		size_t node = nodes[curInd];
		int state = traverser.visit_up(tree, node);
		switch(state){
			case TRAVERSE_CHILDREN:{
				size_t numChildren = tree.num_child_nodes(node);
				const size_t* child = tree.child_node_ids(node);
				for(size_t i = 0; i < numChildren; ++i)
					nodes.push_back(child[i]);
			}break;
			case ABORT_TRAVERSAL:
				return;
			default: break;
		}

		++curInd;
	}

//	traverse downwards
	curInd = (int)nodes.size() - 1;
	while(curInd >= 0){
		size_t node = nodes[curInd];
		traverser.visit_down(tree, node);
		--curInd;
	}

	traverser.end_traversal(tree);
}

template <typename tree_t, typename traverser_t>
int
TraverseDepthFirstRecursion(const tree_t& tree, traverser_t& traverser, int curNode)
{
	int state = traverser.visit_up(tree, curNode);
	switch(state){
		case TRAVERSE_CHILDREN:{
			size_t numChildren = tree.num_child_nodes(curNode);
			const size_t* child = tree.child_node_ids(curNode);
			for(size_t i = 0; i < numChildren; ++i){
				int childState = TraverseDepthFirstRecursion(tree, traverser, child[i]);
				if(childState == ABORT_TRAVERSAL){
					traverser.visit_down(tree, curNode);
					return ABORT_TRAVERSAL;
				}
			}
		} break;
	}
	traverser.visit_down(tree, curNode);
	return state;
}

template <typename tree_t, typename traverser_t>
void TraverseDepthFirst(const tree_t& tree, traverser_t& traverser)
{
	traverser.begin_traversal(tree);
	TraverseDepthFirstRecursion(tree, traverser, 0);
	traverser.end_traversal(tree);
}

}// end of namespace

#endif
