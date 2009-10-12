//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m12 d04

#ifndef __H__LIB_GRID__KD_TREE__
#define __H__LIB_GRID__KD_TREE__

#include "lib_grid/grid/grid.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////////////////////////////
//KDTREE

////////////////////////////////////////////////////////////////////////
///	used by KDTreeStatic
class KDVertexDistance
{
	public:
		KDVertexDistance()	{}
		KDVertexDistance(VertexBase* vrt, float nDistSQ) : vertex(vrt), distSQ(nDistSQ)	{}

		VertexBase*		vertex;
		number			distSQ;
};

////////////////////////////////////////////////////////////////////////
///	used by KDTreeStatic
enum KDSplitDimension
{
	KDSD_CIRCULAR,
	KDSD_LARGEST
};

typedef std::list<KDVertexDistance> KDVertexDistanceList;

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
template <class TPositionAttachment, int numDimensions = 3, class TVector = vector3 >
class KDTreeStatic
{
	public:
		typedef std::vector<VertexBase*> VertexVec;

		class Node
		{
			public:
				Node() : m_pvVertices(NULL)	{m_pChild[0] = m_pChild[1] = NULL;}
				~Node()	{clear();}

				void clear();

				Node*		m_pChild[2];	//	0: pos, 1: neg
				float		m_fSplitValue;
				int			m_iSplitDimension;
				VertexVec*	m_pvVertices;
		};

	//	the functions
		KDTreeStatic() : m_pGrid(NULL)	{};

		void clear();

		bool create_from_grid(Grid& grid, VertexBaseIterator vrtsBegin, VertexBaseIterator vrtsEnd,
								TPositionAttachment& aPos, int maxTreeDepth, int splitThreshold,
								KDSplitDimension splitDimension = KDSD_LARGEST);

		bool get_neighbourhood(std::list<VertexBase*>& vrtsOut,
								typename TPositionAttachment::ValueType& pos, int numClosest);

		bool get_points_in_box(std::list<VertexBase*>& vrtsOut,
								const TVector& boxMin, const TVector& boxMax);

		Node* get_root()	{return &m_parentNode;}
		
		void get_leafs(std::vector<Node*>& vLeafsOut);
		
	protected:
		bool get_points_in_box(std::list<VertexBase*>& vrtsOut, Node* pNode,
								const TVector& boxMin, const TVector& boxMax);

		void neighbourhood(KDVertexDistanceList& vrtsOut, Node* pNode, TVector& pos, int numClosest);

		bool create_barycentric(VertexBaseIterator vrts_begin, VertexBaseIterator vrts_end,
								int numVertices, Node* pNode, int actDimension, int maxTreeDepth);

		int get_largest_dimension(VertexBaseIterator vrts_begin, VertexBaseIterator vrts_end);

		int get_next_split_dimension(int actSplitDimension, VertexBaseIterator vrts_begin,
										VertexBaseIterator vrts_end);
		
		void get_leafs_recursive(std::vector<Node*>& vLeafsOut, Node* pNode);

	//	members
		Grid*	m_pGrid;
		Grid::VertexAttachmentAccessor<TPositionAttachment>	m_aaPos;
		int		m_iSplitThreshold;
		Node	m_parentNode;
		KDSplitDimension	m_splitDimension;	//	how is the next split dimension choosen?

	//	some helper vars for neighbourhood search
		int		m_numNeighboursFound;
		float	m_maxDistSQ;
};

}//	end of namespace

////////////////////////////////////////////////////////////////////////
//	include implementation
#include "kd_tree_static_impl.hpp"

#endif
