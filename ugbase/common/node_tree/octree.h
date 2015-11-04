#ifndef __H__UG__NODE_GRAPH__OCTREE__
#define __H__UG__NODE_GRAPH__OCTREE__

#include "common/math/ugmath.h"
#include "collision_tree_root_node.h"
#include "collision_edges_node.h"
#include "collision_triangles_node.h"

namespace ug{
namespace node_tree
{
////////////////////////////////////////////////////////////////////////
//	CreateOctree
///	Creates an Octree from a list of elements
/**
 * Given a list of points and a list of element-indices (describing edges
 * or triangles) that reference those points, this method constructs
 * a CollisionTree in the form of an octree and returns the root-node
 * of the tree.
 *
 * Arrays passed to this method will be copied to the resulting tree
 * and may thus be deleted after the method returned.
 *
 * \param points: an array of points which are referenced by the triangleInds.
 *
 * \param numPoints: The number of points in the points array.
 *
 * \param elemInds: Array of indices which reference the points array.
 *				'numIndsPerElem' consecutive indices form an element.
 *				Note that the size of elemInds has to be numElemInds.
 *
 * \param numElemInds: The number of elements in the elemInds array.
 *
 * \param numIndsPerElem: The number of points required to describe an element.
 *					Please note, that currently only 2 or 3 are supported.
 *
 * \param elemIDs: Array of size numElemInds / numIndsPerElem, that allows
 *				to specify an identifier associated with each element in the oct-tree.
 *				This parameter may be NULL.
 *
 * \param maxDepth: Maximal depth of the generated trees.
 *
 * \param elemThreshold: The minimal number of elements that leads to an
 *					additional refinement of a node.
 *
 * \param bLoose:
 */
SPCollisionTreeRootNode 
CreateOctree(vector3* points, size_t numPoints,
			  int* elemInds, size_t numElemInds, int numIndsPerElem,
			  CollisionElementID* elemIDs, int maxDepth,
			  int elemThreshold, bool bLoose);


////////////////////////////////////////////////////////////////////////
//	CreateSubOctrees
///	Creates sub-trees of an Octree.
/**
 * This method is mostly used internally during oct-tree-construction.
 * Given a list of points and a element-index-list (edges or triangles),
 * it constructs sub-trees of the specified parent-node.
 * This method recursivly calls itself.
 *
 * Make sure that all passed points lie inside the parent-nodes box.
 *
 * Arrays passed to this method will be copied to the resulting tree
 * and may thus be deleted after the method returned.
 *
 *
 * \param parentNode: The tree-node into which created trees shall be inserted.
 *
 * \param points: an array of points which are referenced by the triangleInds.
 *
 * \param numPoints: The number of points in the points array.
 *
 * \param elemInds: Array of indices which reference the points array.
 *				'numIndsPerElem' consecutive indices form an element.
 *				Note that the size of elemInds has to be numElemInds.
 *
 * \param numElemInds: The number of elements in the elemInds array.
 *
 * \param numIndsPerElem: The number of points required to describe an element.
 *					Please note, that currently only 2 and 3 are supported.
 *
 * \param elemIDs: Array of size numElemInds / numIndsPerElem, that allows
 *				to specify an identifier associated with each element in the oct-tree.
 *				This parameter may be NULL.
 *
 * \param maxDepth: Maximal depth of the generated trees.
 *
 * \param elemThreshold: The minimal number of elements that leads to an
 *					additional refinement of a node.
 *
 * \param bLoose: if true, elems will be assigned to the box that contains
 *				their center. Each element will be assigned to exactly one
 *				leaf-node in this case. Note that the bounding-boxes of a
 *				loose tree normally overlap each other.
 *
 *				if false, elems will be assigned to each leaf-node, whose
 *				box intersects the element. Each element may be assigned to
 *				multiple leaf-nodes in this case. However the boxes of
 *				leaf-nodes will not overlap in this case.
 */
void CreateSubOctrees(BoxedGroupNode* parentNode,
					  vector3* points, size_t numPoints,
			  		  int* elemInds, size_t numElemInds, int numIndsPerElem,
			  		  CollisionElementID* elemIDs, int maxDepth,
					  int elemThreshold, bool bLoose);
								   
}//	end of namesapce node_tree
}//	end of namespace ug

#endif
