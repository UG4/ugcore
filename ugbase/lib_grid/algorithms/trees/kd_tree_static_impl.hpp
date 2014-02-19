//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m12 d04

#ifndef __H__LIB_GRID__KD_TREE_IMPL__
#define __H__LIB_GRID__KD_TREE_IMPL__

#include <list>
#include <vector>
#include "lib_grid/grid/grid.h"
#include "common/math/ugmath.h"

namespace ug
{

typedef std::list<KDVertexDistance> KDVertexDistanceList;

template<class TPositionAttachment, int numDimensions, class TVector>
void
KDTreeStatic<TPositionAttachment, numDimensions, TVector>
::Node
::clear()
{
	if(m_pChild[0])
		delete m_pChild[0];
	if(m_pChild[1])
		delete m_pChild[1];
	m_pChild[0] = m_pChild[1] = NULL;
	if(m_pvVertices)
		delete m_pvVertices;
	m_pvVertices = NULL;
}

template<class TPositionAttachment, int numDimensions, class TVector>
void
KDTreeStatic<TPositionAttachment, numDimensions, TVector>::
clear()
{
	m_pGrid = NULL;
	m_parentNode.clear();
}

template<class TPositionAttachment, int numDimensions, class TVector>
template <class TVrtIterator>
bool
KDTreeStatic<TPositionAttachment, numDimensions, TVector>::
create_from_grid(Grid& grid, TVrtIterator vrtsBegin, TVrtIterator vrtsEnd,
					TPositionAttachment& aPos, int maxTreeDepth,
					int splitThreshold, KDSplitDimension splitDimension)
{
	clear();
	m_pGrid = &grid;
	m_aaPos.access(grid, aPos);
	m_iSplitThreshold = splitThreshold;
	m_splitDimension = splitDimension;	//	how the split dimensions are chosen
	return create_barycentric(vrtsBegin, vrtsEnd, grid.num_vertices(), &m_parentNode, 0, maxTreeDepth);
}

template<class TPositionAttachment, int numDimensions, class TVector>
bool
KDTreeStatic<TPositionAttachment, numDimensions, TVector>::
get_neighbourhood(std::list<Vertex*>& vrtsOut, typename TPositionAttachment::ValueType& pos, int numClosest)
{
	vrtsOut.clear();
	m_numNeighboursFound = 0;
	KDVertexDistanceList pdList;

	neighbourhood(pdList, &m_parentNode, pos, numClosest);
	for(KDVertexDistanceList::iterator iter = pdList.begin(); iter != pdList.end(); iter++)
		vrtsOut.push_back((*iter).vertex);
	return true;
}

template<class TPositionAttachment, int numDimensions, class TVector>
bool
KDTreeStatic<TPositionAttachment, numDimensions, TVector>::
get_points_in_box(std::list<Vertex*>& vrtsOut, const TVector& boxMin, const TVector& boxMax)
{
	vrtsOut.clear();
	return get_vertices_in_box(vrtsOut, &m_parentNode, boxMin, boxMax);
}

template<class TPositionAttachment, int numDimensions, class TVector>
void
KDTreeStatic<TPositionAttachment, numDimensions, TVector>::
get_leafs(std::vector<Node*>& vLeafsOut)
{
	vLeafsOut.clear();
	get_leafs_recursive(vLeafsOut, &m_parentNode);
}


////////////////////////////////////////////////////////////////////////
//	protected
template<class TPositionAttachment, int numDimensions, class TVector>
bool
KDTreeStatic<TPositionAttachment, numDimensions, TVector>::
get_points_in_box(std::list<Vertex*>& vrtsOut, Node* pNode, const TVector& boxMin, const TVector& boxMax)
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

template<class TPositionAttachment, int numDimensions, class TVector>
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
	if(pNode->m_pChild[0] || pNode->m_pChild[1])	//either both or none are NULL
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

template<class TPositionAttachment, int numDimensions, class TVector>
template <class TVertexIterator>
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

template<class TPositionAttachment, int numDimensions, class TVector>
template <class TVertexIterator>
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

template<class TPositionAttachment, int numDimensions, class TVector>
template <class TVertexIterator>
int
KDTreeStatic<TPositionAttachment, numDimensions, TVector>::
get_next_split_dimension(int actSplitDimension, TVertexIterator vrts_begin, TVertexIterator vrts_end)
{
	switch(m_splitDimension)
	{
		case KDSD_LARGEST:
			return get_largest_dimension(vrts_begin, vrts_end);
			break;
		case KDSD_CIRCULAR:
			return (actSplitDimension+1) % numDimensions;
			break;

	}
//	default: SD_CIRCULAR
	return (actSplitDimension+1) % numDimensions;
}

template<class TPositionAttachment, int numDimensions, class TVector>
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
			if(pNode->m_pChild[i] != NULL)
				get_leafs_recursive(vLeafsOut, pNode->m_pChild[i]);
		}
	}
}

}//	end of namespace

#endif
