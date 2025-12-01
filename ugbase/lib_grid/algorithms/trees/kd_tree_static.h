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

#ifndef __H__LIB_GRID__KD_TREE__
#define __H__LIB_GRID__KD_TREE__

#include "lib_grid/grid/grid.h"

namespace ug
{

///	\addtogroup lib_grid_algorithms_trees
///	@{

////////////////////////////////////////////////////////////////////////////////////////////////
//KDTREE

////////////////////////////////////////////////////////////////////////
///	used by KDTreeStatic
class KDVertexDistance
{
	public:
		KDVertexDistance() = default;
		KDVertexDistance(Vertex* vrt, float nDistSQ) : vertex(vrt), distSQ(nDistSQ)	{}

		Vertex* vertex;
		number distSQ;
};

////////////////////////////////////////////////////////////////////////
///	used by KDTreeStatic
enum class KDSplitDimension
{
	KDSD_CIRCULAR,
	KDSD_LARGEST
};

using KDVertexDistanceList = std::list<KDVertexDistance>;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	KDTreeStatic
///	organizes vertices in a binary-tree structure. Only for static use!
/**
 * A kd-tree allows you to find vertices close to a given position in O(log(n)).
 *
 * This kd-tree should be only used for static geometry. If you intend to
 * add or delete vertices after creation, KDTreeStatic is not suited for your
 * needs.
 *
 * This class should be replaced by a dynamic kd-tree, which is capable of
 * dynamic auto-balancing.
 */
template <typename TPositionAttachment, int numDimensions = 3, typename TVector = vector3 >
class KDTreeStatic
{
	public:
		using VertexVec = std::vector<Vertex*>;

		class Node
		{
			public:
				Node() : m_pvVertices(nullptr)	{m_pChild[0] = m_pChild[1] = nullptr;}
				~Node()	{clear();}

				void clear();

				Node*		m_pChild[2];	//	0: pos, 1: neg
				float		m_fSplitValue;
				int			m_iSplitDimension;
				VertexVec*	m_pvVertices;
		};

	//	the functions
		KDTreeStatic() : m_pGrid(nullptr)	{};

		void clear();

		template <typename TVrtIterator>
		bool create_from_grid(Grid& grid, TVrtIterator vrtsBegin, TVrtIterator vrtsEnd,
								TPositionAttachment& aPos, int maxTreeDepth, int splitThreshold,
								KDSplitDimension splitDimension = KDSplitDimension::KDSD_LARGEST);

		template <typename TVrtIterator>
		bool create_from_grid(Grid& grid, TVrtIterator vrtsBegin, TVrtIterator vrtsEnd,
								Grid::VertexAttachmentAccessor<TPositionAttachment> aaPos,
								int maxTreeDepth, int splitThreshold,
								KDSplitDimension splitDimension = KDSplitDimension::KDSD_LARGEST);

		bool get_neighbourhood(std::vector<Vertex*>& vrtsOut,
								typename TPositionAttachment::ValueType& pos, int numClosest);

		bool get_points_in_box(std::vector<Vertex*>& vrtsOut,
								const TVector& boxMin, const TVector& boxMax);

		Node* get_root() {return &m_parentNode;}
		
		void get_leafs(std::vector<Node*>& vLeafsOut);
		
	protected:
		bool get_points_in_box(std::vector<Vertex*>& vrtsOut, Node* pNode,
								const TVector& boxMin, const TVector& boxMax);

		void neighbourhood(KDVertexDistanceList& vrtsOut, Node* pNode, TVector& pos, int numClosest);

		template <typename TVertexIterator>
		bool create_barycentric(TVertexIterator vrts_begin, TVertexIterator vrts_end,
								int numVertices, Node* pNode, int actDimension, int maxTreeDepth);

		template <typename TVertexIterator>
		int get_largest_dimension(TVertexIterator vrts_begin, TVertexIterator vrts_end);

		template <typename TVertexIterator>
		int get_next_split_dimension(int actSplitDimension, TVertexIterator vrts_begin,
										TVertexIterator vrts_end);
		
		void get_leafs_recursive(std::vector<Node*>& vLeafsOut, Node* pNode);

	//	members
		Grid* m_pGrid;
		Grid::VertexAttachmentAccessor<TPositionAttachment>	m_aaPos;
		int m_iSplitThreshold;
		Node m_parentNode;
		KDSplitDimension m_splitDimension;	//	how is the next split dimension choosen?

	//	some helper vars for neighbourhood search
		int m_numNeighboursFound;
		float m_maxDistSQ;
};

/// @}

}//	end of namespace

////////////////////////////////////////////////////////////////////////
//	include implementation
#include "kd_tree_static_impl.hpp"

#endif
