/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_GRID__KD_TREE_IMPL__
#define __H__LIB_GRID__KD_TREE_IMPL__

#include <list>
#include <vector>

#include "lib_grid/grid/grid.h"
// #include "common/math/ugmath.h"

namespace ug {

using KDVertexDistanceList = std::list<KDVertexDistance>;

template <typename TPositionAttachment, int numDimensions, typename TVector>
void
KDTreeStatic<TPositionAttachment, numDimensions, TVector>
::Node
::clear()
{
	if(m_pChild[0])
		delete m_pChild[0];
	if(m_pChild[1])
		delete m_pChild[1];
	m_pChild[0] = m_pChild[1] = nullptr;
	if(m_pvVertices)
		delete m_pvVertices;
	m_pvVertices = nullptr;
}

template <typename TPositionAttachment, int numDimensions, typename TVector>
void
KDTreeStatic<TPositionAttachment, numDimensions, TVector>::
clear()
{
	m_pGrid = nullptr;
	m_parentNode.clear();
}

template<typename TPositionAttachment, int numDimensions, typename TVector>
template <typename TVrtIterator>
bool
KDTreeStatic<TPositionAttachment, numDimensions, TVector>::
create_from_grid(Grid& grid, TVrtIterator vrtsBegin, TVrtIterator vrtsEnd,
					TPositionAttachment& aPos, int maxTreeDepth,
					int splitThreshold, KDSplitDimension splitDimension)
{
	Grid::VertexAttachmentAccessor<TPositionAttachment> aaPos(grid, aPos);
	return create_from_grid(grid, vrtsBegin, vrtsEnd, aaPos, maxTreeDepth,
							splitThreshold, splitDimension);
}

template <typename TPositionAttachment, int numDimensions, typename TVector>
template <typename TVrtIterator>
bool
KDTreeStatic<TPositionAttachment, numDimensions, TVector>::
create_from_grid(Grid& grid, TVrtIterator vrtsBegin, TVrtIterator vrtsEnd,
				 Grid::VertexAttachmentAccessor<TPositionAttachment> aaPos,
				 int maxTreeDepth, int splitThreshold,
				 KDSplitDimension splitDimension)
{
	clear();
	m_pGrid = &grid;
	m_aaPos = aaPos;
	m_iSplitThreshold = splitThreshold;
	m_splitDimension = splitDimension;	//	how the split dimensions are chosen
	return create_barycentric(vrtsBegin, vrtsEnd, grid.num_vertices(), &m_parentNode, 0, maxTreeDepth);
}

template<typename TPositionAttachment, int numDimensions, typename TVector>
bool
KDTreeStatic<TPositionAttachment, numDimensions, TVector>::
get_neighbourhood(std::vector<Vertex*>& vrtsOut, typename TPositionAttachment::ValueType& pos, int numClosest)
{
	vrtsOut.clear();
	m_numNeighboursFound = 0;
	KDVertexDistanceList pdList;

	neighbourhood(pdList, &m_parentNode, pos, numClosest);
	for(KDVertexDistanceList::iterator iter = pdList.begin(); iter != pdList.end(); iter++)
		vrtsOut.push_back((*iter).vertex);
	return true;
}

template<typename TPositionAttachment, int numDimensions, typename TVector>
bool
KDTreeStatic<TPositionAttachment, numDimensions, TVector>::
get_points_in_box(std::vector<Vertex*>& vrtsOut, const TVector& boxMin, const TVector& boxMax)
{
	vrtsOut.clear();
	return get_vertices_in_box(vrtsOut, &m_parentNode, boxMin, boxMax);
}

template<typename TPositionAttachment, int numDimensions, typename TVector>
void
KDTreeStatic<TPositionAttachment, numDimensions, TVector>::
get_leafs(std::vector<Node*>& vLeafsOut)
{
	vLeafsOut.clear();
	get_leafs_recursive(vLeafsOut, &m_parentNode);
}


////////////////////////////////////////////////////////////////////////
//	protected
template<typename TPositionAttachment, int numDimensions, typename TVector>
bool
KDTreeStatic<TPositionAttachment, numDimensions, TVector>::
get_points_in_box(std::vector<Vertex*>& vrtsOut, Node* pNode, const TVector& boxMin, const TVector& boxMax)
{
//	check if we have reached a leaf
	if(pNode->m_pvVertices)
	{
	//	we're in a leaf -> add the points that are inseide the box to pointsOut
		for(VertexVec::iterator iter = pNode->m_pvVertices->begin(); iter != pNode->m_pvVertices->end(); iter++)
		{
			bool bAdd = true;
			for(int i = 0; i < numDimensions; ++i)
			{
				if((m_aaPos[*iter].coord(i) < boxMin[i]) || (m_aaPos[*iter].coord(i) > boxMax[i]))
				{
					bAdd = false;
					break;
				}
			}
			if(bAdd)
				vrtsOut.push_back(*iter);
		}
	//	done
		return true;
	}


//	we haven't reached a leaf yet
//	check in wich subnodes we have to descent
	if(pNode->m_pChild[0] || pNode->m_pChild[1])
	{
		bool bSuccess = true;

		if(boxMin[pNode->m_iSplitDimension] < pNode->m_fSplitValue)
			bSuccess &= get_points_in_box(vrtsOut, pNode->m_pChild[1], boxMin, boxMax);
		if(boxMax[pNode->m_iSplitDimension] >= pNode->m_fSplitValue)
			bSuccess &= get_points_in_box(vrtsOut, pNode->m_pChild[0], boxMin, boxMax);

		return bSuccess;
	}
	return false;
}

template<typename TPositionAttachment, int numDimensions, typename TVector>
void
KDTreeStatic<TPositionAttachment, numDimensions, TVector>::
neighbourhood(KDVertexDistanceList& vrtsOut, Node* pNode, TVector& pos, int numClosest)
{
//	if there are points in this node add them to the list (if they match certain criterias)
	if(pNode->m_pvVertices)
	{
	//	loop through the points associated with the current node
		for(VertexVec::iterator iter = pNode->m_pvVertices->begin(); iter != pNode->m_pvVertices->end(); iter++)
		{
		//	calculate the distance of the current point and the test point
			float distSQ = 0;
			for(int i = 0; i < numDimensions; ++i)
			{
				float t = pos[i] - m_aaPos[*iter].coord(i);
				distSQ += (t*t);
			}

			int bInsert = 0;
			if(m_numNeighboursFound == numClosest)
			{
				if(distSQ < m_maxDistSQ)
					bInsert = 2;
			}
			else
				bInsert = 1;

			if(bInsert)
			{
				KDVertexDistanceList::iterator pdIter;
			//	first we will insert the point into pointlist
				for(pdIter = vrtsOut.begin(); pdIter != vrtsOut.end(); pdIter++)
				{
					if(distSQ < (*pdIter).distSQ)
					{
						pdIter = vrtsOut.insert(pdIter, KDVertexDistance(*iter, distSQ));
						m_numNeighboursFound++;
						break;
					}
				}

				if(pdIter == vrtsOut.end())
				{
					vrtsOut.push_back(KDVertexDistance(*iter, distSQ));
					m_numNeighboursFound++;
				}
				m_maxDistSQ = vrtsOut.back().distSQ;

				if(m_numNeighboursFound > numClosest)
				{
				//	erase the first one
					KDVertexDistanceList::iterator tmpIter = vrtsOut.end();
					tmpIter--;
					vrtsOut.erase(tmpIter);
					m_numNeighboursFound--;
				//	find the new maxDistSQ
					m_maxDistSQ = vrtsOut.back().distSQ;

				}
			}
		}
	}
//	if the node has children visit them
	if(pNode->m_pChild[0] || pNode->m_pChild[1])	//either both or none are nullptr
	{
	//	check in wich subnode the specified point lies.
		int bestNodeIndex;
		if(pos[pNode->m_iSplitDimension] >= pNode->m_fSplitValue)
			bestNodeIndex = 0;
		else
			bestNodeIndex = 1;

	//	first we will visit the better ranked node (just a guess)
		if(pNode->m_pChild[bestNodeIndex])
			neighbourhood(vrtsOut, pNode->m_pChild[bestNodeIndex], pos, numClosest);

	//	if we found less than numClosest we got to continue searching anyway
	//	else if we may find points on the other side of the splitting plane that are closer than the points we found till now we'll have to visit the worse ranked node too.
		if(pNode->m_pChild[1 - bestNodeIndex])
		{
			if(m_numNeighboursFound < numClosest)
				neighbourhood(vrtsOut, pNode->m_pChild[1 - bestNodeIndex], pos, numClosest);
			else
			{
				float t = pos[pNode->m_iSplitDimension] - pNode->m_fSplitValue;
				if(t*t < m_maxDistSQ)
					neighbourhood(vrtsOut, pNode->m_pChild[1 - bestNodeIndex], pos, numClosest);
			}
		}
	}
}

template<typename TPositionAttachment, int numDimensions, typename TVector>
template <typename TVertexIterator>
bool
KDTreeStatic<TPositionAttachment, numDimensions, TVector>::
create_barycentric(TVertexIterator vrts_begin, TVertexIterator vrts_end, int numVertices,
					Node* pNode, int actDimension, int maxTreeDepth)
{
//	check if we are in a leaf
	{
		if((maxTreeDepth < 1) || (numVertices <= m_iSplitThreshold))
		{
		//	we are
			pNode->m_pvVertices = new VertexVec;
		//	loop through the points and add them to the node
			for(TVertexIterator iter = vrts_begin; iter != vrts_end; iter++)
				pNode->m_pvVertices->push_back(*iter);
		//	we're done
			return true;
		}
	}

//	loop through the points and calculate the barycentre
	float barycentre = 0;
	{
		for(TVertexIterator iter = vrts_begin; iter != vrts_end; iter++)
			barycentre += m_aaPos[*iter].coord(actDimension);
		barycentre /= (float)numVertices;
	}

//	fill the lists for the poitive and negative subnodes
	std::list<Vertex*> lstPos;
	std::list<Vertex*> lstNeg;
	int numPos = 0;
	int numNeg = 0;
	{
		for(TVertexIterator iter = vrts_begin; iter != vrts_end; iter++, numPos++, numNeg++)
		{
			if(m_aaPos[*iter].coord(actDimension) >= barycentre)
				lstPos.push_back(*iter);
			else
				lstNeg.push_back(*iter);
		}
	}
//	create the subnodes
	bool bSuccess = true;
	{
		pNode->m_iSplitDimension = actDimension;
		pNode->m_fSplitValue = barycentre;

		for(int i = 0; i < 2; ++i)
			pNode->m_pChild[i] = new Node;
	//	the positive one
		if(!lstPos.empty())
			bSuccess &= create_barycentric(lstPos.begin(), lstPos.end(), numPos, pNode->m_pChild[0], get_next_split_dimension(actDimension, lstPos.begin(), lstPos.end()), maxTreeDepth - 1);
	//	the negative one
		if(!lstNeg.empty())
			bSuccess &= create_barycentric(lstNeg.begin(), lstNeg.end(), numNeg, pNode->m_pChild[1], get_next_split_dimension(actDimension, lstNeg.begin(), lstNeg.end()), maxTreeDepth - 1);
	}
//	done
	return bSuccess;
}

template<typename TPositionAttachment, int numDimensions, typename TVector>
template <typename TVertexIterator>
int
KDTreeStatic<TPositionAttachment, numDimensions, TVector>::
get_largest_dimension(TVertexIterator vrts_begin, TVertexIterator vrts_end)
{
	using namespace std;
	TVector boxMin;
	TVector boxMax;

	TVertexIterator iter = vrts_begin;
//	assign initial values
	{
		for(int i = 0; i < numDimensions; i++)
			boxMin[i] = boxMax[i] = m_aaPos[*iter].coord(i);
	}
	iter++;

//	loop through the points and calculate the bounding box
	for(; iter != vrts_end; iter++)
	{
		for(int i = 0; i < numDimensions; ++i)
		{
			boxMin[i] = min(boxMin[i], m_aaPos[*iter].coord(i));
			boxMax[i] = max(boxMax[i], m_aaPos[*iter].coord(i));
		}
	}
//	calculate extension in each dimension
	TVector extension;
	{
		for(int i = 0; i < numDimensions; ++i)
			extension[i] = boxMax[i] - boxMin[i];
	}
//	return the index of the largest extension
	int bCI = 0;
	{
		for(uint i = 1; i < TVector::Size; ++i)
		{
			if(extension[i] > extension[bCI])
				bCI = (int)i;
		}
	}
	return bCI;
}

template<typename TPositionAttachment, int numDimensions, typename TVector>
template <typename TVertexIterator>
int
KDTreeStatic<TPositionAttachment, numDimensions, TVector>::
get_next_split_dimension(int actSplitDimension, TVertexIterator vrts_begin, TVertexIterator vrts_end)
{
	switch(m_splitDimension)
	{
		case KDSplitDimension::KDSD_LARGEST:
			return get_largest_dimension(vrts_begin, vrts_end);
			break;
		case KDSplitDimension::KDSD_CIRCULAR:
			return (actSplitDimension+1) % numDimensions;
			break;

	}
//	default: SD_CIRCULAR
	return (actSplitDimension+1) % numDimensions;
}

template<typename TPositionAttachment, int numDimensions, typename TVector>
void
KDTreeStatic<TPositionAttachment, numDimensions, TVector>::
get_leafs_recursive(std::vector<Node*>& vLeafsOut, Node* pNode)
{
//	check whether the node is a leaf
	if(!(pNode->m_pChild[0] || pNode->m_pChild[1]))
	{
	//	it is a leaf. push it to the list
		vLeafsOut.push_back(pNode);
	}
	else
	{
		for(int i = 0; i < 2; ++i)
		{
			if(pNode->m_pChild[i] != nullptr)
				get_leafs_recursive(vLeafsOut, pNode->m_pChild[i]);
		}
	}
}

}//	end of namespace

#endif
