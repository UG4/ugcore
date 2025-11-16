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

#ifndef __H__UG__ntree_impl__
#define __H__UG__ntree_impl__

#include <cassert>
#include "ntree.h"

namespace ug{

template <int tree_dim, int world_dim, class elem_t, class common_data_t>
ntree<tree_dim, world_dim, elem_t, common_data_t>::
ntree() :
	m_warningsEnabled (true)
{
	m_nodes.resize(1);
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
void ntree<tree_dim, world_dim, elem_t, common_data_t>::
clear_nodes()
{
	m_nodes.clear();
	m_nodes.resize(1);
	m_nodes[0].level = 0;
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
void ntree<tree_dim, world_dim, elem_t, common_data_t>::
clear()
{
	clear_nodes();
	m_entries.clear();
	m_numDelayedElements = 0;
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
void ntree<tree_dim, world_dim, elem_t, common_data_t>::
set_desc(const NTreeDesc& desc)
{
	m_desc = desc;
}

template <int tree_dim, int world_dim, class elem_t, class common_data_t>
const NTreeDesc& ntree<tree_dim, world_dim, elem_t, common_data_t>::
desc() const
{
	return m_desc;
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
	clear_nodes();
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

	if(m_nodes[nodeIndex].numEntries <= 1)
		return;

	if(m_nodes[nodeIndex].level >= m_desc.maxDepth){
		if(m_warningsEnabled){
			UG_LOG("WARNING in ntree::split_leaf_node(): maximum tree depth "
				<< m_desc.maxDepth << " reached. No further splits are performed for "
				" this node. Note that too many elements per node may lead to performance issues.\n"
				<< "  Number of elements in this node: " << m_nodes[nodeIndex].numEntries << std::endl
				<< "  Corner coordinates of this node: " << m_nodes[nodeIndex].tightBox << std::endl);
		}
		return;
	}

	if(m_nodes[nodeIndex].childNodeInd[0] != s_invalidIndex)
		return;

	const size_t firstChild = m_nodes.size();
	m_nodes.resize(firstChild + s_numChildren);

//	ATTENTION: Be careful not to resize m_nodes while using node, since this would invalidate the reference!
	Node& node = m_nodes[nodeIndex];

//	calculate center of mass and use the traits class to split the box of
//	the current node into 's_numChildren' child boxes. Each child box thereby
//	spanned by one of the corners of the original box and 'centerOfMass'.
	vector_t centerOfMass = calculate_center_of_mass(node);
	box_t childBoxes[s_numChildren];
	traits::split_box(childBoxes, node.tightBox, centerOfMass);

//	iterate over all entries in the current node and assign them to child nodes.
	size_t numEntriesAssigned = 0;
	for(size_t entryInd = node.firstEntryInd; entryInd != s_invalidIndex;){
		Entry& entry = m_entries[entryInd];
		size_t nextEntryInd = entry.nextEntryInd;

		size_t i_child;
		vector_t center;
		traits::calculate_center(center, entry.elem, m_commonData);
		for(i_child = 0; i_child < s_numChildren; ++i_child){
			if(traits::box_contains_point(childBoxes[i_child], center)){
				add_entry_to_node(m_nodes[firstChild + i_child], entryInd);
				++numEntriesAssigned;
				break;
			}
		}
		/*-- For debugging only: --*
		if(i_child == s_numChildren){
			UG_LOG ("ERROR in ntree::split_leaf_node(): Element with center @ " << center
				<< " does not belong to any child of the box " << node.tightBox << std::endl);
		}
		 *--*/

		entryInd = nextEntryInd;
	}

//	all elements of the current box now should be assigned to child boxes.
//	we thus clear element lists and entry-count from the current node.
	UG_COND_THROW(numEntriesAssigned != node.numEntries, "Couldn't find a matching "
				  "child node for some elements during split_leaf_node in "
				  "ntree::split_leaf_node");

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
	
//todo: Solve box-growing in a more portable way!
	typename traits::vector_t offset;
	offset = SMALL;
	traits::grow_box(node.looseBox, node.looseBox, offset);
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
		return elem_iterator_t(nullptr, s_invalidIndex);
	return elem_iterator_t(&m_entries.front(), m_nodes[nodeId].firstEntryInd);
}


template <int tree_dim, int world_dim, class elem_t, class common_data_t>
typename ntree<tree_dim, world_dim, elem_t, common_data_t>::elem_iterator_t
ntree<tree_dim, world_dim, elem_t, common_data_t>::
elems_end(size_t nodeId) const
{
	return elem_iterator_t(nullptr, s_invalidIndex);
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
