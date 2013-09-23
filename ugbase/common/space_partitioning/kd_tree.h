// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Sep 3, 2013


///	JUST AN EARLY MOCKUP - no implementation is yet provided...

#ifndef __H__UG__kd_tree__
#define __H__UG__kd_tree__

namespace ug{

enum KDTreeSplitStrategy
{
	KDTSS_CIRCULAR,
	KDTSS_LARGEST
};

struct KDTreeDesc{
	KDTreeDesc() : dim(3), splitThreshold(16), strategy(KDTSS_LARGEST)	{}
	KDTreeDesc(int nDim, int nSplitThreshold, KDTreeSplitStrategy nStrategy) :
		dim(nDim),
		splitThreshold(nSplitThreshold),
		strategy(nStrategy)	{}

	int dim;
	int splitThreshold;
	KDTreeSplitStrategy	strategy;
};


///	A KDTree is a space-partitioning tree which
template <class point_t, class data_t, class real_t = float>
class KDTree
{
	public:
	///	An Element represents a point with an associated data value.
		class Element{
			friend class KDTree;
			public:
				point_t		point;
				data_t		data;

			private:
			/** index into m_elements. -1: no next element.
			 * Used to create a linked list for each node.*/
				int	nextElemInd;
		};


	///	sets the balancing-parameters of the tree. Triggers a rebalance if a value changed.
		void set_desc(const KDTreeDesc& desc);

	///	returns true if the tree is empty
		bool empty() const;

	///	returns the number of elements in the tree
		size_t size() const;

	///	returns the number of elements which are scheduled for delayed insertion
	/**	\note delayed elements are inserted on a call to rebalance.*/
		size_t num_delayed_elements() const;

	///	adds an element to the tree.
	/**	If 'delayInsertion' is disabled (disabled by default), elements are inserted into
	 * their corresponding nodes on the fly. If the node's splitThreshold is surpassed,
	 * the elements of that node will be sorted into new child nodes automatically.
	 * Note that this may result in an unbalanced tree. You may call rebalance to
	 * reorder elements to receive a balanced tree.
	 * Alternatively you may add_element with 'delayInsertion = true'. All those elements
	 * will be collected and won't be inserted into the tree until rebalance is called.*/
		void add_element(const point_t& point, const data_t& data, bool delayInsertion = false);

	///	finds the closest tree-element to a given point in space
	/** \note If the tree is empty, an std::runtime_error will be thrown.
	 * \note The method call*/
		Element closest_element(const point_t& point) const;

	///	returns the n closest tree-elements to a given point in space
	/** \note when the method is done, elemsOut will have the size min(KDTree::size(), n).*/
		void closest_elements(vector<Element>& elemsOut, const point_t& point,
							  size_t n) const;

	///	estimates the balance quality of the tree in the range [0, 1].
	/**	The returned value represents the ratio to the smallest leaf
	 * (lowest level and smallest number of elements) to the maximum number of
	 * elements, that a node on the same level and its children have.*/
		real_t estimate_balance_quality();

	///	rebalances the whole tree
		void rebalance();

	private:
	/**	The tree is built as a hierarchy of nodes.
	 * Leaf nodes contain elements.*/
		struct Node{
			int			childNodeInds[2]; /// < index into m_nodes. -1: no child node.
			int			firstElemInd; ///< index into m_elements. -1: no element
			int			lastElemInd; ///< index into m_elements. -1: no element
			int			numElements; ///< number of elements in the node
			float		splitValue;
			int			splitDim;
		};

	///	splits a node into two child nodes and assigns elements to those children.
	/**	If the node-threshold of a child node is surpassed, then the child will
	 * be splitted recursively.
	 * Make sure that the given node is a leaf node and thus hasn't got children.*/
		void split_leaf_node(Node& node);

	///	returns an index to the leaf-node which contains the given point
		int find_leaf_node(const point_t& point);


		KDTreeDesc				m_desc;
	/// m_nodes[0] is always considered to be the root node.
		std::vector<Node>		m_nodes;
		std::vector<Element>	m_elements;
		size_t					m_numDelayedElements;
};

}// end of namespace

#endif
