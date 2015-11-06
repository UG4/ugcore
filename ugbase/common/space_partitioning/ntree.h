// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Sep 3, 2013

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
	typedef int	real_t;
	typedef int	vector_t;
	typedef int	box_t;

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

	NTreeDesc() : maxDepth(24), splitThreshold(16)	{}
};


///	The n-tree class can be used to construct e.g. loose quadtrees (tree_dim = 2)
///	or octrees (tree_dim = 3)
/**	world_dim has to be at least as large as tree_dim*/
template <int tree_dim, int world_dim, class TElem, class TCommonData>
class ntree
{
	private:
		struct Entry;

	public:
		typedef TElem										elem_t;
		typedef TCommonData									common_data_t;
		typedef ntree_traits<tree_dim, world_dim, elem_t, common_data_t> traits;
		typedef typename traits::real_t						real_t;
		typedef typename traits::vector_t					vector_t;
		typedef typename traits::box_t						box_t;
		typedef const_ntree_element_iterator<elem_t, Entry>	elem_iterator_t;

		ntree();

		void clear();

	///	sets the balancing-parameters of the tree. Triggers a rebalance
	/** \todo	Only trigger a rebalance if new values require it.*/
		void set_desc(const NTreeDesc& desc);

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
		struct pow<n, 0>	{static const size_t val = 1;};

		template <size_t n, size_t exponent>
		struct pow	{static const size_t val = n * pow<n, exponent - 1>::val;};

	///	the number of children each non-leaf node has
		static const size_t s_numChildren = pow<2, tree_dim>::val;

	///	marks an index as invalid
		static const size_t s_invalidIndex = -1;


	///	An Entry combines an element with the index to the next entry in a leaf-node's entry list.
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

	///	splits a node into two child nodes and assigns entries to those children.
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
};


}// end of namespace


////////////////////////////////////////
//	include implementation
#include "ntree_impl.hpp"

#endif
