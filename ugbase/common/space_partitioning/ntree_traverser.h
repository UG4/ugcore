// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Sep 4, 2013

#ifndef __H__UG__ntree_traverser__
#define __H__UG__ntree_traverser__

namespace ug{
//
//template <class ntree_t>
//class NTreeTraverser_EstimateQuality
//{
//	public:
//		typedef typename ntree_t::traits		traits;
//		typedef typename traits::real_t			real_t;
//		typedef typename traits::vector_t		vector_t;
//		typedef typename traits::box_t			box_t;
//		typedef typename traits::elem_t			elem_t;
//		typedef typename traits::common_data_t	common_data_t;
//
//		NTreeTraverser_EstimateQuality() :
//			m_tree(NULL),
//			m_minLeafNodeLvl(-1),
//			m_curLvl(0),
//			m_minMaxSet(false),
//			m_minNumElems(0),
//			m_maxNumElems(0),
//			m_elemCounter(0),
//			m_searchingMinLeafNode(true)
//		{}
//
//		void traverse_tree(const ntree_t& tree, size_t rootNodeId)
//		{
//			m_searchingMinLeafNode = true;
//			m_minLeafNodeLvl = -1;
//			m_elemCounter = m_minNumElems = m_maxNumElems = 0;
//			m_curLvl = 0;
//			m_minMaxSet = false;
//			m_tree = &tree;
//
//		//	find smallest leaf node
//			m_tree->traverse_node(*this, rootNodeId);
//
//		//	if no smallest leaf-node was found, we'll exit right away...
//			assert(m_minLeafNodeLvl >= 0);
//			if(m_minLeafNodeLvl == -1)
//				return;
//
//		//	now find min- and max-num-elems on m_minLeafNodeLvl (including descendants)
//			m_searchingMinLeafNode = false;
//			m_tree->traverse_node(*this, rootNodeId);
//		}
//
//		void traverse_node(size_t nodeId, const box_t& box,
//						   size_t numElements, bool hasChildren)
//		{
//			if(m_searchingMinLeafNode){
//			//todo:	a breadth first search would be better here...
//				if(!hasChildren){
//					if(m_minLeafNodeLvl == -1)
//						m_minLeafNodeLvl = m_curLvl;
//					else
//						m_minLeafNodeLvl = std::min(m_minLeafNodeLvl, m_curLvl);
//				}
//				else{
//					++m_curLvl;
//					if((m_minLeafNodeLvl == -1) || (m_curLvl < m_minLeafNodeLvl))
//						m_tree->traverse_children(*this, nodeId);
//					--m_curLvl;
//				}
//			}
//			else{
//				if(m_curLvl < m_minLeafNodeLvl){
//					assert(hasChildren);
//					++m_curLvl;
//					m_tree->traverse_children(*this, nodeId);
//					--m_curLvl;
//				}
//				else if(m_curLvl == m_minLeafNodeLvl){
//				//	count children in this node and in all its descendants
//					m_elemCounter = numElements;
//					if(hasChildren){
//						++m_curLvl;
//						m_tree->traverse_children(*this, nodeId);
//						--m_curLvl;
//					}
//					if(m_minMaxSet){
//						m_minNumElems = std::min(m_minNumElems, m_elemCounter);
//						m_maxNumElems = std::max(m_maxNumElems, m_elemCounter);
//					}
//					else{
//						m_minNumElems = m_maxNumElems = m_elemCounter;
//						m_minMaxSet = true;
//					}
//				}
//				else{
//					m_elemCounter += numElements;
//					if(hasChildren){
//						++m_curLvl;
//						m_tree->traverse_children(*this, nodeId);
//						--m_curLvl;
//					}
//				}
//			}
//		}
//
//	///	returns false since no more elements have to be examined
//	/**	This method shouldn't ever be called in this traverser anyways!*/
//		bool traverse_element(const elem_t& elem, const common_data_t& commonData)
//		{return false;}
//
//
//		size_t min_leaf_level() const	{return m_minLeafNodeLvl;}
//		real_t quality() const
//		{
//			if(m_minNumElems == 0)
//				return 0;
//			return (real_t)m_minNumElems / (real_t)m_maxNumElems;
//		}
//
//	private:
//		const ntree_t*	m_tree;
//		int 		m_minLeafNodeLvl;
//		int			m_curLvl;
//		bool		m_minMaxSet;
//		size_t		m_minNumElems;
//		size_t		m_maxNumElems;
//		size_t		m_elemCounter;
//		bool		m_searchingMinLeafNode;
//};

}// end of namespace

#endif
