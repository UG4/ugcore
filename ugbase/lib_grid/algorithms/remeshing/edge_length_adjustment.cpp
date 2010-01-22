//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m01 d22

#include <queue>
#include "lib_grid/lib_grid.h"
#include "edge_length_adjustment.h"

using namespace std;

namespace ug
{
/*
static void AssignFixedVertices(Grid& grid, SubsetHandler& shMarks)
{
//	mark all vertices that lie on a fixed edge as fixed vertex
	if(shMarks.num_subsets() <= RM_FIXED)
		return;

	for(EdgeBaseIterator iter = shMarks.begin<EdgeBase>(RM_FIXED);
		iter != shMarks.end<EdgeBase>(RM_FIXED); ++iter)
	{
		EdgeBase* e = *iter;
		shMarks.assign_subset(e->vertex(0), RM_FIXED);
		shMarks.assign_subset(e->vertex(1), RM_FIXED);
	}
}

static void AssignCreaseVertices(Grid& grid, SubsetHandler& shMarks)
{
//	mark all vertices that lie on a crease and which are not fixed
//	as crease vertices.
	if(shMarks.num_subsets() <= RM_CREASE)
		return;

	for(EdgeBaseIterator iter = shMarks.begin<EdgeBase>(RM_CREASE);
		iter != shMarks.end<EdgeBase>(RM_CREASE); ++iter)
	{
		EdgeBase* e = *iter;
		for(uint i = 0; i < 2; ++i)
			if(shMarks.get_subset_index(e->vertex(i)) != RM_FIXED)
				shMarks.assign_subset(e->vertex(i), RM_CREASE);
	}
}




bool AdjustEdgeLength(Grid& grid, SubsetHandler& shMarks,
					  number minEdgeLen, number maxEdgeLen,
					  int numIterations)
{
//	make sure that faces create associated edges
	if(!grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
	{
		LOG("INFO: auto-enabling FACEOPT_AUTOGENERATE_EDGES in AdjustEdgeLength.\n");
		grid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}

//	replace this by a template parameter
	APosition aPos = aPosition;

//	position attachment
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);

//	check squares
	minEdgeLen *= minEdgeLen;
	maxEdgeLen *= maxEdgeLen;

//	assign vertex marks
	//AssignFixedVertices(grid, shMarks);
	//AssignCreaseVertices(grid, shMarks);

//	we'll store all edges that are candidates in a candidate array
	queue<EdgeBase*> queCandidates;

//	start the main iteration
	for(int iteration = 0; iteration < numIterations; ++iteration)
	{
	//	collect split-candidates
		for(EdgeBaseIterator iter = grid.begin<EdgeBase>();
			iter != grid.end<EdgeBase>(); ++iter)
		{
			EdgeBase* e = *iter;
			if(shMarks.get_subset_index(e) != RM_FIXED)
			{
				if(VecDistanceSq(aaPos[e->vertex(0)], aaPos[e->vertex(1)]) > maxEdgeLen)
					queCandidates.push(e);
			}
		}
	
	//	perform the splits
		//bool performingSplits = !queCandidates.empty();
		while(!queCandidates.empty())
		{
		//	get an candidate
			EdgeBase* e = queCandidates.front();
			queCandidates.pop();

		//	get the center of the edges
			vector3 vCenter = CalculateCenter(e, aaPos);

		//	split the edge
			Vertex* vrt = SplitEdge<Vertex>(grid, e);

		//	assign the new position
			aaPos[vrt] = vCenter;

		//	check if the new edges (edges connected to vCenter)
		//	are candidates again.
			EdgeBaseIterator iterEnd = grid.associated_edges_end(vrt);
			for(EdgeBaseIterator iter = grid.associated_edges_begin(vrt);
				iter != iterEnd; ++iter)
			{
				if(VecDistanceSq(aaPos[e->vertex(0)], aaPos[e->vertex(1)]) > maxEdgeLen)
					queCandidates.push(e);
			}
		}
	}
}
*/
}//	end of namespace
