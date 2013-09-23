// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Sep 4, 2013

#ifndef __H__UG__ntree_traverser__
#define __H__UG__ntree_traverser__

#include <algorithm>
#include <utility>
#include "ntree_traversal.h"

namespace ug{

template <class tree_t>
class Traverser_FindLowestLeafNodeLevel
{
	public:
		Traverser_FindLowestLeafNodeLevel() :
			m_lowestLeafNodeLvl(0)	{}

		void begin_traversal(const tree_t& tree)
		{
			m_lowestLeafNodeLvl = 0;
		}

		bool visit_up(const tree_t& tree, size_t node)
		{
			if(tree.num_child_nodes(node) == 0){
				m_lowestLeafNodeLvl = tree.level(node);
				return ABORT_TRAVERSAL;
			}
			return TRAVERSE_CHILDREN;
		}

		void visit_down(const tree_t&, size_t)	{}

		void end_traversal(const tree_t&)	{}

		size_t result() const {return m_lowestLeafNodeLvl;}

	private:
		size_t m_lowestLeafNodeLvl;
};

template <class tree_t>
size_t FindLowestLeafNodeLevel(const tree_t& tree)
{
	Traverser_FindLowestLeafNodeLevel<tree_t> t;
	TraverseBreadthFirst(tree, t);
	return t.result();
}


///	returns the minimum and maximum number of elements in all subtrees of nodes of the given level
template <class tree_t>
class Traverser_MinMaxNumElements
{
	public:
		Traverser_MinMaxNumElements(size_t lvl) :
			m_lvl(lvl), m_minNumElements(0), m_maxNumElements(0),
			m_elemCount(0), m_firstEval(true)	{}

		void begin_traversal(const tree_t& tree)
		{
			m_minNumElements = m_maxNumElements = 0;
			m_elemCount = 0;
			m_firstEval = true;
		}

		bool visit_up(const tree_t& tree, size_t node)
		{
			if(tree.level(node) == m_lvl)
				m_elemCount = 0;

			if(tree.level(node) >= m_lvl)
				m_elemCount += tree.num_elements(node);

			return TRAVERSE_CHILDREN;
		}

		void visit_down(const tree_t& tree, size_t node)
		{
			if(tree.level(node) == m_lvl){
				if(m_firstEval){
					m_minNumElements = m_maxNumElements= m_elemCount;
					m_firstEval = false;
				}
				else{
					m_minNumElements = std::min(m_minNumElements, m_elemCount);
					m_maxNumElements = std::max(m_maxNumElements, m_elemCount);
				}
			}
		}

		void end_traversal(const tree_t&)	{}

		size_t min_num_elements() const {return m_minNumElements;}
		size_t max_num_elements() const {return m_maxNumElements;}

	private:
		size_t m_lvl;
		size_t m_minNumElements;
		size_t m_maxNumElements;
		size_t m_elemCount;
		bool m_firstEval;
};

template <class tree_t>
std::pair<size_t, size_t> GetMinMaxNumElements(const tree_t& tree, size_t lvl)
{
	Traverser_MinMaxNumElements<tree_t> t(lvl);
	TraverseDepthFirst(tree, t);
	return std::make_pair(t.min_num_elements(), t.max_num_elements());
}

}// end of namespace

#endif
