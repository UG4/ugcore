#ifndef __H__UG__NODE_TREE__COLLSISION_TRIANGLES_NODE__
#define __H__UG__NODE_TREE__COLLSISION_TRIANGLES_NODE__

#include <vector>
#include "node.h"
#include "collision_element_info.h"

namespace ug{
namespace node_tree
{
class CollisionTrianglesNode;

////////////////////////////////////////////////////////////////////////
///	the smartpointer used to encapsulate the node
typedef SmartPtr<CollisionTrianglesNode> SPCollisionTrianglesNode;

////////////////////////////////////////////////////////////////////////
//	CollisionTrianglesNode
///	holds index tuples defining triangles.
/**
 * The index-tuples refer to the next CollisionTreeRootNode, which is
 * higher in the hierarchy.
 *
 * An identifier can be stored with each triangle - either an
 * integer-value or a void-pointer.
 * Normally appears as a subordinate of CollisionTreeRootNode.
 */
class CollisionTrianglesNode : public Node
{
	public:
		static SPCollisionTrianglesNode create();

		virtual ~CollisionTrianglesNode();

		virtual void add_triangle(int ind1, int ind2, int ind3);
		virtual void add_triangle(int ind1, int ind2, int ind3,
								  CollisionElementID triID);
		
	/// pIndices has to be of size numTris*3
		virtual void add_triangles(int* pIndices, size_t numTris);
		
	///	pIndices and pTriIDs have to be of size numTris*3
		virtual void add_triangles(int* pIndices,
								   CollisionElementID* pTriIDs,
								   size_t numTris);
		
		virtual size_t num_triangles() const;
		
		virtual void get_triangle(size_t index, int& ind1Out,
								  int& ind2Out, int& ind3Out) const;

		virtual const int* get_triangles() const;

		virtual void set_triangle_id(size_t triInd,
									 CollisionElementID triID);
		
	/// if no identifier has been set for an edge -1 is returned.
		virtual CollisionElementID get_triangle_id(size_t triInd);

	protected:
		CollisionTrianglesNode();

	protected:
		typedef std::vector<int>				IndexVec;
		typedef std::vector<CollisionElementID>	IDVec;

		IndexVec	m_vTris;
		IDVec		m_vTriIDs;
		
		bool m_bTriangleIDsSupplied;
};

}//	end of namespace node_tree
}//	end of namespace ug

#endif
