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

#ifndef __H__UG__ntree__
#define __H__UG__ntree__

#include "ntree_iterator.h"

namespace ug{

/**	The following methods have to be provided for the given vector-type:
 * \code
 * void VecSet(vector_t& vOut, real_t value);	// sets all components of 'v' to 'value'.
 * void VecAdd(vector_t& vOut, const vector_t& v1, const vector_t& v2); // performs vOut = v1 + v2.
 * void VecScale(vector_t& vOut, const vector_t& v, real_t s); // performs vOut = s * v.
 * \endcode
 */
template <int tree_dim, int world_dim, class elem_t, class common_data_t>
struct ntree_traits
{
	using real_t = int;
	using vector_t = int;
	using box_t = int;

	static void calculate_center(vector_t& centerOut, const elem_t& e,
						  	  	 const common_data_t& commonData);

	static void calculate_bounding_box(box_t& boxOut, const elem_t& e,
									   const common_data_t& commonData);

///	adds the given offset to box.max and subtracts it from box.min
	static void grow_box(box_t& boxOut, const box_t& box,
						 const vector_t& offset);

	static vector_t box_diagonal(const box_t& box);

	static bool box_contains_point(const box_t& box, const vector_t& point);

///	returns true if the given boxes intersect
	static bool box_box_intersection(const box_t& box1, const box_t& box2);

///	returns the smallest box that contains both box1 and box2
	static void merge_boxes(box_t& boxOut, const box_t& box1, const box_t& box2);

///	splits the given box into (2^tree_dim sub-boxes).
/**	The split should be performed so that the given split-point is the only common
 * point of all boxes.
 * The union of all boxes in boxesOut has to overlap the given input-box.*/
	static void split_box(box_t* boxesOut, const box_t& box, const vector_t& splitPoint);

////	todo: the following methods should go into special traverser traits.
///** required for ContainsPointTraverser.*/
	static bool contains_point(const elem_t& e, const vector_t& point,
							   const common_data_t& commonData);
//
///**	required for ClosestEntriesTraverser.*/
//	real_t distance_point_to_entry(const elem_t& e, const vector_t& point,
//									 const elem_data_t& elemData,
//									 const common_data_t& commonData);
//
///**	required for IntersectingEntriesTraverser.*/
//	bool ray_intersects_entry(real_t distOut, const elem_t& elem,
//								const vector_t& rayFrom, const vector_t& rayTo,
//								const elem_data_t& elemData,
//								const common_data_t& commonData);
};




struct NTreeDesc{
	size_t maxDepth;
	size_t splitThreshold;
	//todo:	split-strategy (center-of-space, center-of-mass)

	NTreeDesc() : maxDepth(32), splitThreshold(16)	{}
};


///	The n-tree class can be used to construct space partitioning trees of dimensions 1, 2, and 3.
/** The ntree class provides spatial partitioning through binary- (tree_dim=1),
 * quad- (tree_dim=2), and octrees (tree_dim=3). It is designed as a template class to
 * support both different types of underlying math libraries and also arbitrary objects
 * that shall be partitioned. The element type is specified through TElem. It
 * does not have to fulfill any concepts but should be lightweight, since it is
 * copied during tree creation.
 *
 * The tree furthermore features a common_data object whose type is determined
 * through a template argument again. This common_data object is passed to
 * traversers and traits. Required objects such as data-accessors can be stored
 * in this common_data object as required by a specialization of ntree_traits.
 * It does not have to fulfill any concepts in regard to the ntree itself.
 *
 * Usage:
 *	- Add elements to your tree-instance through the 'add_element' method.
 *	- Call 'rebalance' once all elements have been added.
 *	- Use traversers to query the tree for elements, e.g.
 *	  'Traverser_FindContainingElement', 'Traverser_FindElementsInIntersectingNodes',
 *	  or 'Traverser_RayElementIntersection'.
 *
 * \note	Whether a node is refined into child-nodes during 'rebalance' is
 *			determined by the 'splitThreshold' of the associated 'NTreeDesc'
 *			instance. It can be adjusted through 'set_desc'.
 *			Only if a node contains more than 'splitThreshold' elements, a split
 *			of that node will be performed. Default is 16.
 *
 * \note	Each tree has a maximum tree depth. It defaults to 32 and can be
 *			adjusted through the 'set_desc' method. If this tree depth is
 *			reached, the tree won't be splitted any further and a warning
 *			is printed. Note that a tree of this depth often results from a bad
 *			underlying geometry and may lead to performance issues.		
 *
 * \param tree_dim		Dimension of the tree: binary-trees (tree_dim=1),
 * 						quad-trees (tree_dim=2), and octrees (tree_dim=3)
 *						are supported.
 *
 * \param world_dim		Dimension of the space in which the tree is embedded.
 *						This primarily affects underlying mathematical structures
 *						(i.e. vectors). Please note that world_dim has to be at
 *						least as large as tree_dim.
 *
 * \param TElem			The element type for which the tree is constructed. It
 *						should be lightweight, since it is copied during tree
 *						creation. There are no special concepts that TElem has
 *						to fullfill.
 *
 * \param TCommonData	User provided data that is stored in the tree (one instance only)
 *						but not used by the tree. It is passed to functions in
 *						ntree_traits. The concrete type for TCommonData depends on
 *						the concrete ntree_traits that shall be used.
 */
template <int tree_dim, int world_dim, class TElem, class TCommonData>
class ntree
{
	private:
		struct Entry;

	public:
		using elem_t = TElem;
		using common_data_t = TCommonData;
		using traits = ntree_traits<tree_dim, world_dim, elem_t, common_data_t>;
		using real_t = typename traits::real_t;
		using vector_t = typename traits::vector_t;
		using box_t = typename traits::box_t;
		using elem_iterator_t = const_ntree_element_iterator<elem_t, Entry>;

		ntree();

		void clear();

	///	enabled or disable warning messages
	/** warnings are enabled by default. If a problem is detected and warnings
	 * are enabled, a warning message will be written to stdout.*/
		void enable_warnings(bool enable)	{m_warningsEnabled = enable;}
		bool warnings_enabled () const		{return m_warningsEnabled;}

	///	sets the balancing-parameters of the tree.
	/** \note	The method has no effect until the next call to 'rebalance',
	 *			i.e., it will not affect an existing tree if one has already
	 *			been built.*/
		void set_desc(const NTreeDesc& desc);

	///	returns the balancing-parameters of the tree.
		const NTreeDesc& desc() const;

	///	sets the common-data which the tree passes on to callback methods
		void set_common_data(const common_data_t& commonData);

	///	returns the common-data stored in the tree
		const common_data_t& common_data() const;

	///	returns true if the tree is empty
		bool empty() const;

	///	returns the number of entries in the tree (delayed entries excluded)
		size_t size() const;

	///	returns the number of elements which have been added but are not yet accessible in the tree.
	/**	\note delayed elements are inserted on a call to rebalance.*/
		size_t num_delayed_elements() const;

	///	adds an element to the tree.
	/**	The element will only be scheduled for insertion but won't be inserted
	 * until rebalance is called.
	 * \todo	Allow optional on-the-fly insertion of elements.*/
		void add_element(const elem_t& elem);

	///	rebalances the whole tree
	/**	The method returns false if some error occurred during rebalancing.*/
		void rebalance();

	///	returns the total number of nodes in the tree
		size_t num_nodes() const;

	///	returns the number of children of a node
		size_t num_child_nodes(size_t nodeId) const;

	///	returns an array of child-id's for the given node
	/**	Use ntree::num_children on the node to retrieve the length of the returned array.*/
		const size_t* child_node_ids(size_t nodeId) const;

	///	returns an iterator to the first element of a given node
		elem_iterator_t elems_begin(size_t nodeId) const;

	///	returns an iterator to the end of the element-sequence of a given node
	/**	this iterator points one element behind the last element of the sequence.*/
		elem_iterator_t elems_end(size_t nodeId) const;

	///	returns the number of elements that the given node contains
		size_t num_elements(size_t nodeId) const;

	///	returns the number tree-level in which the node is located
		size_t level(size_t nodeId) const;

	///	returns the smallest box which contains all elements of the given node
		const box_t& bounding_box(size_t nodeId) const;

	private:
	
	///	clear only the nodes
		void clear_nodes ();
		
	///	static template implementation to raise n to the power exponent
		template <size_t n, size_t exponent>
		struct pow;

		template <size_t n>
		struct pow<n, 0>	{static constexpr size_t val = 1;};

		template <size_t n, size_t exponent>
		struct pow	{static constexpr size_t val = n * pow<n, exponent - 1>::val;};

	///	the number of children each non-leaf node has
		static constexpr size_t s_numChildren = pow<2, tree_dim>::val;

	///	marks an index as invalid
		static constexpr size_t s_invalidIndex = -1;


	///	An Entry combines an element with the index to the next entry in a leaf-node's entry list.
	/** Note that exactly one 'Entry' object per element exists. Since 'Entry'
	 * also serves as a linked list, this means that an element can only be contained
	 * in one node at a time.*/
		struct Entry{
				elem_t		elem;
			/** index into m_entries. s_invalidIndex: no next entry.
			 * Used to create a linked list for each leaf-node.*/
				size_t	nextEntryInd;

				Entry(const elem_t& e) :
					elem(e), nextEntryInd(s_invalidIndex)	{}
		};


	/**	The tree is built as a hierarchy of nodes.
	 * Leaf nodes contain entries.*/
		struct Node{
			size_t		childNodeInd[s_numChildren]; /// < index into m_nodes. s_invalidIndex: no child node.
			size_t		firstEntryInd; ///< index into m_entries. s_invalidIndex: no entry
			size_t		lastEntryInd; ///< index into m_entries. s_invalidIndex: no entry
			size_t		numEntries; ///< number of entries in the node
			size_t		level;
			box_t		tightBox; ///< tight bounding box - disjunct partition of the root box
			box_t		looseBox; ///< loose bounding box - contains all bounding boxes of its entries

			Node() : firstEntryInd(s_invalidIndex), lastEntryInd(s_invalidIndex),
					 numEntries(0), level(s_invalidIndex)
			{
				tightBox = looseBox;
				for(size_t i = 0; i < s_numChildren; ++i)
					childNodeInd[i] = s_invalidIndex;
			}

			Node(const Node& n) : firstEntryInd(n.firstEntryInd), lastEntryInd(n.lastEntryInd),
								  numEntries(n.numEntries), level(n.level),
								  tightBox(n.tightBox), looseBox(n.looseBox)
			{
				for(size_t i = 0; i < s_numChildren; ++i)
					childNodeInd[i] = n.childNodeInd[i];
			}
		};

	///	splits a node into 2^tree_dim child nodes and assigns entries to those children.
	/**	If the node-threshold of a child node is surpassed, then the child will
	 * be splitted recursively.
	 * Make sure that the given node is a leaf node and thus hasn't got children.*/
		void split_leaf_node(size_t nodeIndex);

	///	returns an index to the leaf-node which contains the given point
	/**	returns s_invalidIndex if no matching node was found.
	 * Checks the point against the tight bounding-box of each node.*/
		size_t find_leaf_node(const vector_t& point, size_t curNode = 0);

	///	adds an entry to the given node
		void add_entry_to_node(Node& node, size_t entryInd);

	///	updates the loose bounding box of the given node
		void update_loose_bounding_box(Node& node);

	///	calculates the center of mass of a given node
		vector_t calculate_center_of_mass(Node& node);

		NTreeDesc				m_desc;
		common_data_t			m_commonData;
		std::vector<Node>		m_nodes; ///< m_nodes[0] is always considered to be the root node.
		std::vector<Entry>		m_entries;
		size_t					m_numDelayedElements;
		bool					m_warningsEnabled;
};


}// end of namespace


////////////////////////////////////////
//	include implementation
#include "ntree_impl.hpp"

#endif
