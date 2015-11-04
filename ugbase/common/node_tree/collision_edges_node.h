#ifndef __H__UG__NODE_TREE__COLLSISION_EDGES_NODE__
#define __H__UG__NODE_TREE__COLLSISION_EDGES_NODE__

#include <vector>
#include "node.h"
#include "collision_element_info.h"

namespace ug{
namespace node_tree
{
class CollisionEdgesNode;

////////////////////////////////////////////////////////////////////////
///	the smartpointer used to encapsulate the node
typedef SmartPtr<CollisionEdgesNode> SPCollisionEdgesNode;

////////////////////////////////////////////////////////////////////////
//	CollisionEdgesNode
///	holds index pairs defining edges.
/**
 * An identifier can be stored with each edge.
 * Normally appears as a subordinate of CollisionTreeRootNode.
 */
class CollisionEdgesNode : public Node
{
	public:
		static SPCollisionEdgesNode create();

		virtual ~CollisionEdgesNode();

		virtual void add_edge(int ind1, int ind2);
		virtual void add_edge(int ind1, int ind2,
							  CollisionElementID edgeID);
		
	/// pIndices has to be of size numEdges*2
		virtual void add_edges(int* pIndices, int numEdges);
		
	///	pIndices has to be of size numEdges*2
		virtual void add_edges(int* pIndices,
							   CollisionElementID* pEdgeIDs,
							   int numEdges);
							   
		virtual int num_edges() const;
		
		virtual void get_edge(int index, int& ind1Out, int& ind2Out) const;
		
		virtual const int* get_edges() const;

		virtual void set_edge_id(int edgeIndex,
								 CollisionElementID edgeID);
		
	/// if no identifier has been set for an edge -1 is returned.
		virtual CollisionElementID get_edge_id(int edgeIndex);

	protected:
		CollisionEdgesNode();

	protected:
		typedef std::vector<int>				IndexVec;
		typedef std::vector<CollisionElementID>	IDVec;

		IndexVec	m_vEdges;
		IDVec		m_vEdgeIDs;

		bool m_bEdgeIDsSupplied;
};

}//	end of namespace node_tree
}//	end of namespace ug

#endif
