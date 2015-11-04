#ifndef __H__UG__ntree_traversal__
#define __H__UG__ntree_traversal__

namespace ug{

enum TraversalStates{
	DONT_TRAVERSE_CHILDREN,
	TRAVERSE_CHILDREN,
	ABORT_TRAVERSAL
};

template <class tree_t, class traverser_t>
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

template <class tree_t, class traverser_t>
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

template <class tree_t, class traverser_t>
void TraverseDepthFirst(const tree_t& tree, traverser_t& traverser)
{
	traverser.begin_traversal(tree);
	TraverseDepthFirstRecursion(tree, traverser, 0);
	traverser.end_traversal(tree);
}

}// end of namespace

#endif
