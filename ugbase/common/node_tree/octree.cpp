/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#include <cassert>
#include <vector>
#include "octree.h"

using namespace std;

namespace ug{
namespace node_tree
{

////////////////////////////////////////////////////////////////////////
/**	calculates the bounding box around a set of points*/
static void
CalculateBoundingBox(vector3& boxMinOut, vector3& boxMaxOut,
					const vector3* points, size_t numPoints)
{
	if(numPoints < 1){
		boxMinOut = boxMaxOut = vector3(0, 0, 0);
		return;
	}
	
	boxMinOut = boxMaxOut = points[0];
	
	for(size_t i = 0; i < numPoints; ++i){
		for(size_t j = 0; j < 3; ++j)
		{
			const number& coord = points[i][j];
			if(coord < boxMinOut[j])
				boxMinOut[j]  = coord;
			else if(coord > boxMaxOut[j])
				boxMaxOut[j]  = coord;
		}
	}
}

////////////////////////////////////////////////////////////////////////
/**	calculates the bounding box around a set of elements which are
 *	given by a set of indices referencing a set of points. This method
 *	is particularily useful if the elements only reference a subset
 * 	of the given points. If not the bounding box should be calculated
 *	directly on the set of points.
 */
static void
CalculateBoundingBox(vector3& boxMinOut, vector3& boxMaxOut,
					const vector3* points,
					const int* inds, size_t numInds)
{
	if(numInds < 1){
		boxMinOut = boxMaxOut = vector3(0, 0, 0);
		return;
	}
	
	boxMinOut = boxMaxOut = points[inds[0]];
	
	for(size_t i = 0; i < numInds; ++i){
		for(size_t j = 0; j < 3; ++j)
		{
			const number& coord = points[inds[i]][j];
			if(coord < boxMinOut[j])
				boxMinOut[j]  = coord;
			else if(coord > boxMaxOut[j])
				boxMaxOut[j]  = coord;
		}
	}
}

////////////////////////////////////////////////////////////////////////
/**	Increases the size of the bounding box by scaling around its center*/
static void GrowBox(vector3& minInOut, vector3& maxInOut, number scaleFac)
{
	vector3 halfDiag;
	VecSubtract(halfDiag, maxInOut, minInOut);
	VecScale(halfDiag, halfDiag, 0.5);
	VecScale(halfDiag, halfDiag, scaleFac);
	
	vector3 center;
	VecAdd(center, minInOut, maxInOut);
	VecScale(center, center, 0.5);
	
	VecAdd(maxInOut, center, halfDiag);
	VecSubtract(minInOut, center, halfDiag);
}

////////////////////////////////////////////////////////////////////////
//	CreateOctree
SPCollisionTreeRootNode 
CreateOctree(vector3* points, size_t numPoints,
			  int* elemInds, size_t numElemInds, int numIndsPerElem,
			  CollisionElementID* elemIDs, int maxDepth,
			  int elemThreshold, bool bLoose)
{
//	if there are no elements, we'll return an empty node
	if(numElemInds == 0)
		return SPCollisionTreeRootNode(nullptr);
		
//	if the number of indices per elements is not right, we'll return, too.
	if(numIndsPerElem < 2 || numIndsPerElem > 3)
		return SPCollisionTreeRootNode(nullptr);

//	calculate the bounding-box of the tree
	vector3 boxMin, boxMax;
	CalculateBoundingBox(boxMin, boxMax, points, numPoints);
	GrowBox(boxMin, boxMax, 1.001);
	
//	create the root-node
	SPCollisionTreeRootNode spRootNode = CollisionTreeRootNode::create();
	
//	add the given points to the root-node
	spRootNode->add_points(points, numPoints);
	
//	set the bounding-box of the root-node
	spRootNode->set_box(boxMin, boxMax);
	
//	create the sub-trees
	CreateSubOctrees(spRootNode.get(), points, numPoints, elemInds,
					numElemInds, numIndsPerElem, elemIDs, maxDepth,
					elemThreshold, bLoose);

//	done. return the root-node
	return spRootNode;
}

////////////////////////////////////////////////////////////////////////
//	CreateSubOctrees
void CreateSubOctrees(BoxedGroupNode* parentNode,
					  vector3* points, size_t numPoints,
			  		  int* elemInds, size_t numElemInds, int numIndsPerElem,
			  		  CollisionElementID* elemIDs, int maxDepth,
					  int elemThreshold, bool bLoose)
{
//	the maximal number of indices per element.
//	has to be adjusted if new element-types are supported.
	assert(numIndsPerElem > 1 && numIndsPerElem <= 3 &&
			"unsupported number of indices per element.");

//	the number of elements
	int numElems = numElemInds / numIndsPerElem;
	
//	check if subnodes have to be created
	if((maxDepth > 0) && (numElems > elemThreshold))
	{
	//	yes, we have to create subnodes.
		vector3 boxMin, boxMax;
		boxMin = parentNode->min_corner();
		boxMax = parentNode->max_corner();
		
		vector3 boxSize;
		boxSize.x() = boxMax.x() - boxMin.x();
		boxSize.y() = boxMax.y() - boxMin.y();
		boxSize.z() = boxMax.z() - boxMin.z();
		
	//	we'll store index-lists of subtrees in those vectors
		vector<int> 				vNewElems;
		vector<CollisionElementID>	vNewIDs;
//TODO:	check whether vNewElems.reserve amd vNewIDs.reserve would result in a speedup.

	//	if a loose tree shall be generated, we'll have to record
	//	whether an edge has already been assigned or not.
		vector<char> elemAssigned;
		if(bLoose)
			elemAssigned.resize(numElems, 0);

	//	iterate through each sub-tree
		for(size_t subNodeInd = 0; subNodeInd < 8; ++subNodeInd)
		{
		//	clear arrays
			vNewElems.clear();
			vNewIDs.clear();

		//	create the box of the subnode
			vector3 subBoxMin, subBoxMax;

		//	compute bounding box from ParentBox and index
			subBoxMin.x() = boxMin.x() + .5 * boxSize.x() * (subNodeInd % 2);
			subBoxMin.y() = boxMin.y() + .5 * boxSize.y() * (int)(subNodeInd / 4);
			subBoxMin.z() = boxMin.z() + .5 * boxSize.z() * (int)((subNodeInd % 4) / 2);
			subBoxMax.x() = subBoxMin.x() + .5 * boxSize.x();
			subBoxMax.y() = subBoxMin.y() + .5 * boxSize.y();
			subBoxMax.z() = subBoxMin.z() + .5 * boxSize.z();
			
		//	find the elements that are inside the box.
		//	different approaches have to be taken, depending on
		//	whether a loose tree shall be created or not.
			if(bLoose){
				GrowBox(subBoxMin, subBoxMax, 1.001);
			//	create a loose tree.
			//	check for each element whether the center lies in the box or not.
				for(size_t i = 0; i < numElemInds; i += numIndsPerElem)
				{
					size_t elemIndex = i / numIndsPerElem;
					
					if(!elemAssigned[elemIndex])
					{
					//	calculate the center
						vector3 center(0, 0, 0);
						for(int j = 0; j < numIndsPerElem; ++j){
							VecAdd(center, center, points[elemInds[i + j]]);
						}
						
						VecScale(center, center, 1. / numIndsPerElem);
					
					//	check whether the center lies in the bounding box
						if(BoxBoundProbe(center, subBoxMin, subBoxMax))
						{
						//	it does. Mark the element as assigned and
						//	add it to the new-lists.
							elemAssigned[elemIndex] = 1;
							
							for(int j = 0; j < numIndsPerElem; ++j)
								vNewElems.push_back(elemInds[i + j]);
								
							if(elemIDs)
								vNewIDs.push_back(elemIDs[elemIndex]);
						}
					}
				}
				
			//	all elements for this sub-tree have been collected.
			//	Since it is a loose tree, we'll have to readjust the box
				CalculateBoundingBox(subBoxMin, subBoxMax, points,
									&vNewElems.front(), vNewElems.size());
									
				GrowBox(subBoxMin, subBoxMax, 1.001);
			}
			else{
//TODO:	check whether moving the element-iteration into the switch would
//		result in a noticeable speedup (I don't think so!).

			//	the tree is not loose. We have to find all elements that
			//	intersect the box.
				for(size_t i = 0; i < numElemInds; i += numIndsPerElem)
				{
				//	check whether elements are lines or triangles
					bool intersects = false;
					switch(numIndsPerElem)
					{
					case 2:	// elements are line-segments
						intersects = LineBoxIntersection(points[elemInds[i]],
														points[elemInds[i + 1]],
														subBoxMin, subBoxMax);
						break;
						
					case 3: // elements are triangles
						intersects = TriangleBoxIntersection(
										points[elemInds[i]], points[elemInds[i + 1]],
										points[elemInds[i + 2]], subBoxMin, subBoxMax);
						break;
					default:
						intersects = false;
						break;
					}
					
				//	if we found an intersection, the element has to be assigned
				//	to the subtree, which is associated with the box.
					if(intersects){
						for(int j = 0; j < numIndsPerElem; ++j)
							vNewElems.push_back(elemInds[i + j]);
						
						if(elemIDs)
							vNewIDs.push_back(elemIDs[i / numIndsPerElem]);
					}
				}
			}

		//	all elements that go into this subtree have been found.
		//	now create the subtree
			SPBoxedGroupNode spSubNode = BoxedGroupNode::create();
			spSubNode->set_box(subBoxMin, subBoxMax);
			parentNode->add_child(spSubNode);

			if(elemIDs)
				CreateSubOctrees(spSubNode.get(), points, numPoints,
								&vNewElems.front(), vNewElems.size(), numIndsPerElem,
								&vNewIDs.front(), maxDepth - 1, elemThreshold, bLoose);
			else
				CreateSubOctrees(spSubNode.get(), points, numPoints,
								&vNewElems.front(), vNewElems.size(), numIndsPerElem,
								nullptr, maxDepth - 1, elemThreshold, bLoose);
		}
	}
	else{
	//	there are less elements thatn elemThreshold or maxDepth is reached.
	//	We will now create the tree-node that holds the elements.
	//	If there are no elements, we can return immediatly.
		if(numElemInds == 0)
			return;
			
	//	We have to distinguish between the different element-types
		switch(numIndsPerElem){
		case 2://	edges
			{
				SPCollisionEdgesNode spNode = CollisionEdgesNode::create();
				if(elemIDs)
					spNode->add_edges(elemInds, elemIDs, numElemInds / 2);
				else
					spNode->add_edges(elemInds, numElemInds / 2);
					
				parentNode->add_child(spNode);
			}
			break;
			
		case 3://	triangles
			{
				SPCollisionTrianglesNode spNode = CollisionTrianglesNode::create();
				if(elemIDs)
					spNode->add_triangles(elemInds, elemIDs, numElemInds / 3);
				else
					spNode->add_triangles(elemInds, numElemInds / 3);

				parentNode->add_child(spNode);
			}
			break;
			
		default:
			break;
		}
	}
}


}//	end of namesapce node_tree
}//	end of namespace ug

