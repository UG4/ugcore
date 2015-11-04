// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Sep 4, 2013 (d,m,y)

#ifndef __H__UG__ntree_impl__
#define __H__UG__ntree_impl__

#include <cassert>
#include "ntree.h"

namespace ug{

template <int tree_dim, int world_dim, class elem_t, class common_data_t>
ntree<tree_dim, world_dim, elem_t, common_data_t>::
ntree()
{
	m_nodes.resize(1);
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
void ntree<tree_dim, world_dim, elem_t, common_data_t>::
clear()
{
	m_nodes.clear();
	m_nodes.resize(1);
	m_nodes[0].level = 0;
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
void ntree<tree_dim, world_dim, elem_t, common_data_t>::
set_desc(const NTreeDesc& desc)
{
	m_desc = desc;
	rebalance();
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
void ntree<tree_dim, world_dim, elem_t, common_data_t>::
set_common_data(const common_data_t& commonData)
{
	m_commonData = commonData;
}

template <int tree_dim, int world_dim, class elem_t, class common_data_t>
const common_data_t& ntree<tree_dim, world_dim, elem_t, common_data_t>::
common_data() const
{
	return m_commonData;
}

template <int tree_dim, int world_dim, class elem_t, class common_data_t>
bool ntree<tree_dim, world_dim, elem_t, common_data_t>::
empty() const
{
	return size() == 0;
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
size_t ntree<tree_dim, world_dim, elem_t, common_data_t>::
size() const
{
	return m_entries.size() - m_numDelayedElements;
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
size_t ntree<tree_dim, world_dim, elem_t, common_data_t>::
num_delayed_elements() const
{
	return m_numDelayedElements;
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
void ntree<tree_dim, world_dim, elem_t, common_data_t>::
add_element(const elem_t& elem)
{
//	size_t entryInd = m_entries.size();
	m_entries.push_back(Entry(elem));

	++m_numDelayedElements;

//todo	If we allowed for immediate on-the-fly insertion of elements, the following
//		code could serve as an implementation base. Note, that it is not yet complete.
//		One would have to update parent boxes if a child box changes etc.
//		It would make sense to also allow for the specification of a
//		minimum root-bounding-box, since else most insertions would trigger a rebalance.
//	{
//		vector_t center;
//		traits::calculate_center(center, elem, m_commonData);
//		size_t nodeInd = find_leaf_node(center);
//		if(nodeInd == s_invalidIndex){
//			++m_numDelayedElements;
//			rebalance();
//		}
//		else{
//			Node& node = m_nodes[nodeInd];
//			add_entry_to_node(node, entryInd);
//			if(node.numEntries >= m_desc.splitThreshold)
//				split_leaf_node(nodeInd);
//			else{
//			//todo	the update could be much more efficient (simply unite the
//			//		old bounding box with the new element's bounding box.
//				update_loose_bounding_box(node);
//			}
//		}
//	}
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
void ntree<tree_dim, world_dim, elem_t, common_data_t>::
rebalance()
{
//	push all elements into the root node, calculate its bounding box
//	and call split_leaf_node if the element threshold is surpassed.
	clear();
	Node& root = m_nodes.back();
	if(!m_entries.empty()){
		root.firstEntryInd = 0;
		root.lastEntryInd = m_entries.size() - 1;
		root.numEntries = m_entries.size();
		for(size_t i = 0; i < m_entries.size(); ++i)
			m_entries[i].nextEntryInd = i+1;
		m_entries.back().nextEntryInd = s_invalidIndex;

	//	the tight bounding box and the loose bounding box of the root node
	//	should be the same.
		update_loose_bounding_box(root);
		root.tightBox = root.looseBox;

		if(root.numEntries >= m_desc.splitThreshold)
			split_leaf_node(0);
	}
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
void ntree<tree_dim, world_dim, elem_t, common_data_t>::
split_leaf_node(size_t nodeIndex)
{

	assert((m_nodes[nodeIndex].numEntries >= 2)
		   && "Not enough elements in a node during split_leaf_node");


	if(m_nodes[nodeIndex].childNodeInd[0] != s_invalidIndex)
		return;

	const size_t firstChild = m_nodes.size();
	m_nodes.resize(firstChild + s_numChildren);

//	ATTENTION: Be careful not to resize m_nodes while using node, since this would invalidate the reference!
	Node& node = m_nodes[nodeIndex];

	vector_t centerOfMass = calculate_center_of_mass(node);
	box_t childBoxes[s_numChildren];
	traits::split_box(childBoxes, node.tightBox, centerOfMass);

	size_t numEntriesAssigned = 0;
	for(size_t entryInd = node.firstEntryInd; entryInd != s_invalidIndex;){
		Entry& entry = m_entries[entryInd];
		size_t nextEntryInd = entry.nextEntryInd;

		vector_t center;
		traits::calculate_center(center, entry.elem, m_commonData);
		for(size_t i_child = 0; i_child < s_numChildren; ++i_child){
			if(traits::box_contains_point(childBoxes[i_child], center)){
				add_entry_to_node(m_nodes[firstChild + i_child], entryInd);
				++numEntriesAssigned;
				break;
			}
		}

		entryInd = nextEntryInd;
	}

	assert((numEntriesAssigned == node.numEntries)
		   && "Couldn't find a matching child node for some elements during split_leaf_node:");

	node.firstEntryInd = node.lastEntryInd = s_invalidIndex;
	node.numEntries = 0;

	for(size_t i_child = 0; i_child < s_numChildren; ++i_child){
		node.childNodeInd[i_child] = firstChild + i_child;
		Node& childNode = m_nodes[firstChild + i_child];
		childNode.level = node.level + 1;
		childNode.tightBox = childBoxes[i_child];
		update_loose_bounding_box(childNode);
	}

//	since split_leaf_node resizes m_nodes and since this invalidates any references
//	to m_nodes, we perform the recursion in a last step
	for(size_t i_child = 0; i_child < s_numChildren; ++i_child){
		size_t childNodeInd = firstChild + i_child;
		if(m_nodes[childNodeInd].numEntries >= m_desc.splitThreshold)
			split_leaf_node(childNodeInd);
	}
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
size_t ntree<tree_dim, world_dim, elem_t, common_data_t>::
find_leaf_node(const vector_t& point, size_t curNode)
{
	Node& n = m_nodes[curNode];
	if(traits::box_contains_point(n.tightBox, point)){
		if(n.childNodeInd[0] != s_invalidIndex){
			for(size_t i = 0; i < s_numChildren; ++i){
				size_t result = find_leaf_node(point, n.childNodeInd[i]);
				if(result != s_invalidIndex)
					return result;
			}
		}
		else{
			return curNode;
		}
	}
	return s_invalidIndex;
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
void ntree<tree_dim, world_dim, elem_t, common_data_t>::
add_entry_to_node(Node& node, size_t entryInd)
{
	if(node.firstEntryInd == s_invalidIndex){
		assert(node.numEntries == 0);
		node.firstEntryInd = node.lastEntryInd = entryInd;
	}
	else{
		m_entries[node.lastEntryInd].nextEntryInd = entryInd;
		node.lastEntryInd = entryInd;
	}

	m_entries[entryInd].nextEntryInd = s_invalidIndex;
	++node.numEntries;
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
void ntree<tree_dim, world_dim, elem_t, common_data_t>::
update_loose_bounding_box(Node& node)
{
	size_t entryInd = node.firstEntryInd;
	if(entryInd == s_invalidIndex){
		node.looseBox = node.tightBox;
		return;
	}

	Entry& entry = m_entries[entryInd];
	traits::calculate_bounding_box(node.looseBox, entry.elem, m_commonData);

	entryInd = entry.nextEntryInd;
	while(entryInd != s_invalidIndex){
		Entry& entry = m_entries[entryInd];
		box_t nbox;
		traits::calculate_bounding_box(nbox, entry.elem, m_commonData);
		traits::merge_boxes(node.looseBox, node.looseBox, nbox);
		entryInd = entry.nextEntryInd;
	}
	
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
typename ntree<tree_dim, world_dim, elem_t, common_data_t>::vector_t
ntree<tree_dim, world_dim, elem_t, common_data_t>::
calculate_center_of_mass(Node& node)
{
	vector_t centerOfMass;
	VecSet(centerOfMass, 0);

	for(size_t entryInd = node.firstEntryInd; entryInd != s_invalidIndex;){
		Entry& entry = m_entries[entryInd];

		vector_t center;
		traits::calculate_center(center, entry.elem, m_commonData);
		VecAdd(centerOfMass, centerOfMass, center);
		entryInd = entry.nextEntryInd;
	}

	if(node.numEntries > 0)
		VecScale(centerOfMass, centerOfMass, (real_t)1 / (real_t)node.numEntries);

	return centerOfMass;
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
size_t ntree<tree_dim, world_dim, elem_t, common_data_t>::
num_nodes() const
{
	return m_nodes.size();
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
size_t ntree<tree_dim, world_dim, elem_t, common_data_t>::
num_child_nodes(size_t nodeId) const
{
	assert(nodeId < m_nodes.size());
	if(m_nodes[nodeId].childNodeInd[0] == s_invalidIndex)
		return 0;
	return s_numChildren;
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
const size_t* ntree<tree_dim, world_dim, elem_t, common_data_t>::
child_node_ids(size_t nodeId) const
{
	assert(nodeId < m_nodes.size());
	return m_nodes[nodeId].childNodeInd;
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
typename ntree<tree_dim, world_dim, elem_t, common_data_t>::elem_iterator_t
ntree<tree_dim, world_dim, elem_t, common_data_t>::
elems_begin(size_t nodeId) const
{
	assert(nodeId < m_nodes.size());
	if(m_entries.empty())
		return elem_iterator_t(NULL, s_invalidIndex);
	return elem_iterator_t(&m_entries.front(), m_nodes[nodeId].firstEntryInd);
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
typename ntree<tree_dim, world_dim, elem_t, common_data_t>::elem_iterator_t
ntree<tree_dim, world_dim, elem_t, common_data_t>::
elems_end(size_t nodeId) const
{
	return elem_iterator_t(NULL, s_invalidIndex);
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
size_t ntree<tree_dim, world_dim, elem_t, common_data_t>::
num_elements(size_t nodeId) const
{
	assert(nodeId < m_nodes.size());
	return m_nodes[nodeId].numEntries;
}

template <int tree_dim, int world_dim, class elem_t, class common_data_t>
size_t ntree<tree_dim, world_dim, elem_t, common_data_t>::
level(size_t nodeId) const
{
	assert(nodeId < m_nodes.size());
	return m_nodes[nodeId].level;
}

template <int tree_dim, int world_dim, class elem_t, class common_data_t>
const typename ntree<tree_dim, world_dim, elem_t, common_data_t>::box_t&
ntree<tree_dim, world_dim, elem_t, common_data_t>::
bounding_box(size_t nodeId) const
{
	assert(nodeId < m_nodes.size());
	return m_nodes[nodeId].looseBox;
}

}// end of namespace

#endif
