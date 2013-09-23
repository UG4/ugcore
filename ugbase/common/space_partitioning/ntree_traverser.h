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

		int visit_up(const tree_t& tree, size_t node)
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
	Traverser_FindLowestLeafNodeLevel<tree_t> trav;
	TraverseBreadthFirst(tree, trav);
	return trav.result();
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

		int visit_up(const tree_t& tree, size_t node)
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
	Traverser_MinMaxNumElements<tree_t> trav(lvl);
	TraverseDepthFirst(tree, trav);
	return std::make_pair(trav.min_num_elements(), trav.max_num_elements());
}


template <class tree_t>
class Traverser_FindContainingElement
{
	public:
		typedef typename tree_t::elem_t		elem_t;
		typedef typename tree_t::vector_t	vector_t;

		Traverser_FindContainingElement(const vector_t& point) :
			m_point(point),
			m_foundElem(false)
		{}

		void begin_traversal(const tree_t& tree)
		{
			m_foundElem = false;
			m_numElemsChecked = 0;
		}

		int visit_up(const tree_t& tree, size_t node)
		{
		//	if the point doesn't lie in the node's bounding box, we don't have
		//	to check it's elements at all.
			if(!tree_t::traits::box_contains_point(tree.bounding_box(node), m_point))
				return DONT_TRAVERSE_CHILDREN;

		//	first check whether the nodes box contains the given point
			if(tree.num_child_nodes(node) == 0){
			//	iterate over all elements. If an element contains the given point,
			//	we're done and we may return.
				for(typename tree_t::elem_iterator_t iter = tree.elems_begin(node);
					iter != tree.elems_end(node); ++iter)
				{
					++m_numElemsChecked;
					if(tree_t::traits::contains_point(*iter, m_point, tree.common_data())){
						m_foundElem = true;
						m_elem = *iter;
						return ABORT_TRAVERSAL;
					}
				}
			}
			return TRAVERSE_CHILDREN;
		}

		void visit_down(const tree_t&, size_t)	{}

		void end_traversal(const tree_t&)	{}

		bool result(elem_t& foundElemOut) const
		{
			if(m_foundElem)
				foundElemOut = m_elem;
			return m_foundElem;
		}

		size_t num_elems_checked() const	{return m_numElemsChecked;}

	private:
		vector_t	m_point;
		elem_t		m_elem;
		size_t		m_numElemsChecked;
		bool		m_foundElem;
};

template <class tree_t>
bool FindContainingElement(typename tree_t::elem_t& elemOut, const tree_t& tree,
						   const typename tree_t::vector_t& point)
{
	Traverser_FindContainingElement<tree_t> trav(point);
	TraverseDepthFirst(tree, trav);
	//UG_LOG("num elems checked for one pick: " << trav.num_elems_checked() << "\n");
	typename tree_t::elem_t tmpElem;
	if(trav.result(tmpElem)){
		elemOut = tmpElem;
		return true;
	}
	return false;
}

}// end of namespace

#endif
