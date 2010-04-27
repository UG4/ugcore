//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m01 d22

#include <queue>
#include "lib_grid/lib_grid.h"
#include "common/profiler/profiler.h"
#include "simple_grid.h"
#include "edge_length_adjustment.h"

using namespace std;

namespace ug
{

////////////////////////////////////////////////////////////////////////
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



////////////////////////////////////////////////////////////////////////
template <class TVertexPositionAccessor>
number CalculateNormalDot(TriangleDescriptor& td1, TriangleDescriptor& td2,
						  TVertexPositionAccessor& aaPos)
{
	vector3 n1;
	CalculateTriangleNormal(n1, aaPos[td1.vertex(0)],
							aaPos[td1.vertex(1)], aaPos[td1.vertex(2)]);
	vector3 n2;
	CalculateTriangleNormal(n2, aaPos[td2.vertex(0)],
							aaPos[td2.vertex(1)], aaPos[td2.vertex(2)]);
	return VecDot(n1, n2);
}


////////////////////////////////////////////////////////////////////////
//	CalculateCurvature
template <class TAAPosVRT>
number CalculateMinCurvature(Grid& grid, SubsetHandler& shMarks,
							VertexBase* vrt, TAAPosVRT& aaPos)
{
//TODO:	check whether static vNormals brings any benefits.
//TODO:	special cases for crease vertices

//	face normals
	static vector<vector3>	vNormals;
//	vertex normal (mean face normal)
	vector3 n(0, 0, 0);
	
//	calculate the normals of associated faces
	vNormals.clear();
	FaceIterator iterEnd = grid.associated_faces_end(vrt);
	for(FaceIterator iter = grid.associated_faces_begin(vrt);
		iter != iterEnd; ++iter)
	{
		vector3 nTmp;
		CalculateNormal(nTmp, *iter, aaPos);
		vNormals.push_back(nTmp);
		VecAdd(n, n, nTmp);
	}

//	the vertex normal
	VecNormalize(n, n);
	
//	get the min dot-product of the vertex normal with associated faces-normals
	number minDot = 1;
	for(size_t i = 0; i < vNormals.size(); ++i)
		minDot = std::min(minDot, VecDot(n, vNormals[i]));

//	done
	return minDot;	
}							

////////////////////////////////////////////////////////////////////////
template <class TAAPosVRT>
number CalculateAverageCurvature(Grid& grid, SubsetHandler& shMarks,
								EdgeBase* e, TAAPosVRT& aaPos)
{
	return 0.5 * (CalculateMinCurvature(grid, shMarks, e->vertex(0), aaPos)
				+ CalculateMinCurvature(grid, shMarks, e->vertex(1), aaPos));
}

////////////////////////////////////////////////////////////////////////
template <class TAAPosVRT>
number CalculateLengthFac(Grid& grid, SubsetHandler& shMarks,
								EdgeBase* e, TAAPosVRT& aaPos)
{
	number lenFac = CalculateAverageCurvature(grid, shMarks, e, aaPos);
	lenFac = (lenFac - 0.95) / 0.05;
	return max(0.25, lenFac);
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
template <class TAAPosVRT, class TAANormVRT, class TAAIntVRT>
bool TrySwap(Grid& grid, EdgeBase* e, TAAPosVRT& aaPos, TAANormVRT& aaNorm,
			TAAIntVRT& aaInt, SubsetHandler* pshMarks = NULL,
			EdgeSelector* pCandidates = NULL)
{
//	swaps are neither allowed for crease edges nor for fixed edges
	if(pshMarks){
		if(pshMarks->get_subset_index(e) == RM_FIXED ||
			pshMarks->get_subset_index(e) == RM_CREASE)
			return false;
	}

//	get the associated faces. we need two of them
	Face* f[2];
	if(GetAssociatedFaces(f, grid, e, 2) != 2)
		return false;

//	make sure both faces are triangles
	if((f[0]->num_vertices() != 3) || (f[1]->num_vertices() != 3))
		return false;

//	create a simple grid
	SimpleGrid sg;
	if(!ObtainSimpleGrid(sg, grid, e->vertex(0), e->vertex(1), 0,
						aaPos, aaNorm, aaInt))
	{
		LOG("ObtainSimpleGrid failed. ignoring edge...\n");
		return false;
	}

//	calculate geometric-approximation-degree and triangle quality
	number approxDeg = GeometricApproximationDegree(sg);
	number shapeDeg = ShapeQualityDegree(sg);

	number smoothDeg = VecDot(sg.triangleNormals[0], sg.triangleNormals[1]);

//	perform a swap on the simple grid
	if(!SwapEdge(sg))
	{
		LOG("swap edge failed...\n");
		return false;
	}

//	calculate new geometric-approximation-degree and triangle quality
	number newApproxDeg = GeometricApproximationDegree(sg);
	number newShapeDeg = ShapeQualityDegree(sg);
	number newSmoothDeg = VecDot(sg.triangleNormals[0], sg.triangleNormals[1]);

//	neither the shapeDeg nor the approxDeg may get too bad.
	if((newApproxDeg < 0.5 * approxDeg) || (newShapeDeg < 0.5 * shapeDeg))
		return false;
	
//	make sure that the swap does not destroy the smoothness
	if(newSmoothDeg < 0.1 * smoothDeg)
		return false;

	//if((newApproxDeg - approxDeg) + 0.5 * (newShapeDeg - shapeDeg) > 0)
	if(newShapeDeg > shapeDeg)
	{
	//	swap the edge
		e = SwapEdge(grid, e);

		if(e){
		//	swap was successful
		//	if pCandidates was specified then add new candidates
			if(pCandidates){
				for(int i = 0; i < 2; ++i)
					pCandidates->select(grid.associated_edges_begin(e->vertex(i)),
										grid.associated_edges_end(e->vertex(i)));
			//	e was selected but is not really a candidate
				pCandidates->deselect(e);
			}
			
			return true;
		}
	}
	
	return false;
}

////////////////////////////////////////////////////////////////////////
template <class TAAPosVRT, class TAANormVRT, class TAAIntVRT>
bool PerformSwaps(Grid& grid, SubsetHandler& shMarks, EdgeSelector& esel,
				TAAPosVRT& aaPos, TAANormVRT& aaNorm, TAAIntVRT& aaInt)
{	
	PROFILE_FUNC();
	LOG("  performing swaps\n");
	int numSwaps = 0;
	
	while(!esel.empty())
	{
		EdgeBase* e = *esel.begin<EdgeBase>();
		esel.deselect(e);
			
		if(TrySwap(grid, e, aaPos, aaNorm, aaInt, &shMarks, &esel))
			++numSwaps;
	}
	LOG("  swaps performed: " << numSwaps << endl);
	
	return true;
}

template <class TAAPosVRT, class TAANormVRT, class TAAIntVRT>
bool TryCollapse(Grid& grid, EdgeBase* e, 
				TAAPosVRT& aaPos, TAANormVRT& aaNorm, 
				TAAIntVRT& aaInt, SubsetHandler* pshMarks = NULL,
				EdgeSelector* pCandidates = NULL)
{
	if(pshMarks)
	{
		SubsetHandler& shMarks = *pshMarks;
	//	collapses are not allowed for fixed edges
		if(shMarks.get_subset_index(e) == RM_FIXED)
			return false;
			
	//	if both endpoints of are fixed vertices then
	//	we may not collapse
		int vrtSI[2];
		vrtSI[0] = shMarks.get_subset_index(e->vertex(0));
		vrtSI[1] = shMarks.get_subset_index(e->vertex(1));
		if((vrtSI[0] == RM_FIXED) && (vrtSI[1] == RM_FIXED))
			return false;

	//	if both endpoints are somehow marked, e has to be a
	//	crease edge
		if((vrtSI[0] != RM_NONE) && (vrtSI[1] != RM_NONE)
			&&	(shMarks.get_subset_index(e) != RM_CREASE))
			return false;
	}

//	check whether the edge can be collapsed
	if(EdgeCollapseIsValid(grid, e))
	{
	//	test the collapse using a simple-grid
		SimpleGrid sg;
		if(!ObtainSimpleGrid(sg, grid, e->vertex(0), e->vertex(1), 1,
							aaPos, aaNorm, aaInt))
		{
			LOG("ObtainSimpleGrid failed. ignoring edge...\n");
			return false;
		}
		
	//	calculate geometric-approximation-degree and triangle quality
		number approxDeg = GeometricApproximationDegree(sg);
		number shapeDeg = ShapeQualityDegree(sg);

	//	perform a swap on the simple grid
		if(!CollapseEdge(sg))
		{
			LOG("collapse edge failed...\n");
			return false;
		}

	//	get the positions of the old endpoints
		static const int numTestPositions = 3;
		int newInd = sg.vertices.size() - 1;
		vector3 v[numTestPositions];
		v[0] = aaPos[e->vertex(0)];
		v[1] = aaPos[e->vertex(1)];
		v[2] = sg.vertices[newInd];
		
	//	we'll compare 3 approximation degrees and three shape degrees
		number newApproxDeg[numTestPositions];
		number newShapeDeg[numTestPositions];
		
	//	check which position is the best
		int bestIndex = -1;
	//	the vertex subset index is used to support marks (crease and fixed vertices)
		int vrtSI[2];
		vrtSI[0] = vrtSI[1] = RM_NONE;

		if(pshMarks)
		{
			vrtSI[0] = pshMarks->get_subset_index(e->vertex(0));
			vrtSI[1] = pshMarks->get_subset_index(e->vertex(1));

			if((vrtSI[0] == RM_FIXED) || ((vrtSI[0] != RM_NONE) && (vrtSI[1] == RM_NONE))){
				bestIndex = 0;
				newApproxDeg[0] = GeometricApproximationDegree(sg);
				newShapeDeg[0] = ShapeQualityDegree(sg);
			}
			else if((vrtSI[1] == RM_FIXED) || ((vrtSI[1] != RM_NONE) && (vrtSI[0] == RM_NONE))){
				bestIndex = 1;
				newApproxDeg[1] = GeometricApproximationDegree(sg);
				newShapeDeg[1] = ShapeQualityDegree(sg);
			}
		}
		
		if(bestIndex == -1){
		//	check all three approximation degrees
			for(int i = 0; i < numTestPositions; ++i){
			//	we'll compute all qualities with the averaged normal
				sg.vertices[newInd] = v[i];
				CalculateTriangleNormals(sg);
				newApproxDeg[i] = GeometricApproximationDegree(sg);
				newShapeDeg[i] = ShapeQualityDegree(sg);
			}
		//	get the best one
			bestIndex = 0;
			for(int i = 1; i < numTestPositions; ++i){
				if(newApproxDeg[i] > newApproxDeg[bestIndex])
					bestIndex = i;
			}
		}
		
	//	if the shape-degree of the collapsed region is too bad, we'll skip the collapse
		if(newShapeDeg[bestIndex] < 0.5 * shapeDeg)
			return false;

	//	if the best approximation degree is not too bad, we'll perform the collapse
		if(newApproxDeg[bestIndex] > 0.8 * approxDeg)
		{						
		//	pick one of the endpoints to be the one that resides
			VertexBase* vrt = e->vertex(1);
			if(vrtSI[0] != RM_NONE)
				vrt = e->vertex(0);
			
		//	collapse the edge
			CollapseEdge(grid, e, vrt);

		//	assign best position
			aaPos[vrt] = v[bestIndex];
		//	assign the normal
			aaNorm[vrt] = sg.vertexNormals[newInd];

			if(pCandidates){
//TODO: all edges that belong to associated faces are new candidates.						
			//	associated edges of vrt are possible new candidates
				pCandidates->select(grid.associated_edges_begin(vrt),
							grid.associated_edges_end(vrt));
			}

			return true;
		}
	}

	return false;
}

////////////////////////////////////////////////////////////////////////
template <class TAAPosVRT, class TAANormVRT, class TAAIntVRT>
bool PerformCollapses(Grid& grid, SubsetHandler& shMarks, EdgeSelector& esel,
					  number minEdgeLen, TAAPosVRT& aaPos, TAANormVRT& aaNorm,
					  TAAIntVRT& aaInt)
{	
	PROFILE_FUNC();
	LOG("  performing collapses\n");
	int numCollapses = 0;
//	compare squares
	minEdgeLen *= minEdgeLen;

	while(!esel.empty())
	{
		EdgeBase* e = *esel.begin<EdgeBase>();
		esel.deselect(e);
		
	//	the higher the curvature the smaller the maxEdgeLen.
	//	minimal lenFac is 0.1
		number lenFac = CalculateLengthFac(grid, shMarks, e, aaPos);

	//	check whether the edge is short enough
		if(VecDistanceSq(aaPos[e->vertex(0)], aaPos[e->vertex(1)]) < lenFac * minEdgeLen)
		{
			if(TryCollapse(grid, e, aaPos, aaNorm, aaInt, &shMarks, &esel))
				++numCollapses;
		}
	}
	LOG("  collapses performed: " << numCollapses << endl);
	
	return true;
}

////////////////////////////////////////////////////////////////////////
template <class TAAPosVRT, class TAANormVRT>
bool TrySplit(Grid& grid, EdgeBase* e, TAAPosVRT& aaPos, TAANormVRT& aaNorm,
			  EdgeSelector* pCandidates = NULL, SubsetHandler* pshMarks = NULL)
{
	bool bCreaseEdge = false;
//	splits are not allowed for fixed edges
	if(pshMarks){
		if(pshMarks->get_subset_index(e) == RM_FIXED)
			return false;
		else if(pshMarks->get_subset_index(e) == RM_CREASE)
			bCreaseEdge = true;
	}

//	get the center of the edges
	vector3 vCenter = CalculateCenter(e, aaPos);

//	the new normal
	vector3 n;
	if(bCreaseEdge){
		VecAdd(n, aaNorm[e->vertex(0)], aaNorm[e->vertex(1)]);
		VecNormalize(n, n);
	}
/*
	int ni = 0;
	if(pshMarks)
		ni = pshMarks->get_subset_index(e);
*/	
//	split the edge
	Vertex* vrt = SplitEdge<Vertex>(grid, e, false);

/*
	if(pshMarks)
		pshMarks->assign_subset(vrt, ni);
*/
//	assign the new position
	aaPos[vrt] = vCenter;

//	assign the new normal. calculate it if required
	if(!bCreaseEdge)
		CalculateVertexNormal(n, grid, vrt, aaPos);			

	aaNorm[vrt] = n;

//	associated edges of vrt are candidates again
	if(pCandidates)
		pCandidates->select(grid.associated_edges_begin(vrt),
							grid.associated_edges_end(vrt));

	return true;
}

////////////////////////////////////////////////////////////////////////
template <class TAAPosVRT, class TAANormVRT>
bool PerformSplits(Grid& grid, SubsetHandler& shMarks, EdgeSelector& esel,
					  number maxEdgeLen, TAAPosVRT& aaPos, TAANormVRT& aaNorm)
{
	PROFILE_FUNC();
//	compare squares
	maxEdgeLen *= maxEdgeLen;

	LOG("  performing splits\n");
	int numSplits = 0;

	while(!esel.empty())
	{
	//	get an edge
		EdgeBase* e = *esel.begin<EdgeBase>();
		esel.deselect(e);
		
	//	the higher the curvature the smaller the maxEdgeLen.
	//	minimal lenFac is 0.1
		number lenFac = CalculateLengthFac(grid, shMarks, e, aaPos);

	//	check whether the edge should be splitted
		if(VecDistanceSq(aaPos[e->vertex(0)], aaPos[e->vertex(1)]) > lenFac * maxEdgeLen)
		{
			if(TrySplit(grid, e, aaPos, aaNorm, &esel, &shMarks))
				++numSplits;
		}
	}

	LOG("  splits performed: " << numSplits << endl);
	return true;
}

////////////////////////////////////////////////////////////////////////
//	relocate point by smoothing
void RelocatePointBySmoothing(vector3& vOut, const vector3&v,
						  std::vector<vector3>& vNodes,
						  size_t numIterations, number stepSize,
						  std::vector<number>* pvWeights = NULL)
{
//TODO:	incorporate weights
//	iterate through all nodes, calculate the direction
//	from v to each node, and add the scaled sum to v.
	vector3 t, vOld;
	vOut = v;
	stepSize /= (number)vNodes.size();

	for(size_t j = 0; j < numIterations; ++j){
		vOld = vOut;
		for(size_t i = 0; i < vNodes.size(); ++i){
			VecSubtract(t, vNodes[i], vOld);
			VecScale(t, t, stepSize);
			VecAdd(vOut, vOut, t);
		}
	}
}



////////////////////////////////////////////////////////////////////////
bool AdjustEdgeLength(Grid& gridOut, SubsetHandler& shOut, SubsetHandler& shMarksOut,
					  Grid& gridIn, SubsetHandler& shIn, SubsetHandler& shMarksIn,
					  number minEdgeLen, number maxEdgeLen, int numIterations,
					  bool projectPoints)
{
	PROFILE_FUNC();
//TODO: check crease-marks and fixed marks.
//		make sure that edge-collapse leaves the grid regular.
//		swaps should be the final operation
//		separate collapses, swaps and splits

//	replace this by a template parameter
	APosition aPos = aPosition;
	
//	we need a reference grid and reference marks
	Grid* pRefGrid = &gridIn;
	SubsetHandler* pRefMarks = &shMarksIn;

//	if the input grid and the output grid are the same, we'll need
//	a temporary reference grid.
//	same goes for marks
	Grid tmpRefGrid;
	SubsetHandler tmpRefMarks;
	if(&gridOut == &gridIn){
		tmpRefGrid = gridIn;
		pRefGrid = &tmpRefGrid;
		
		if(&shMarksOut == &shMarksIn){
			tmpRefMarks.assign_grid(tmpRefGrid);
			tmpRefMarks = shMarksIn;
			pRefMarks = &tmpRefMarks;
		}
	}
	else{
		gridOut = gridIn;
	}

	if(&shMarksOut != &shMarksIn){
		shMarksOut = shMarksIn;
	}

//	initialize shOut
	if(&shOut != &shIn){
		shOut = shIn;
	}

	Grid& grid = gridOut;
	SubsetHandler& shMarks = shMarksOut;
		
//	make sure that grid and pRefGrid have position-attachments
	if(!(grid.has_vertex_attachment(aPos) && pRefGrid->has_vertex_attachment(aPos))){
		UG_LOG("  vertex-position-attachment missing in AdjustEdgeLength. Aborting.\n");
		return false;
	}
	
//	make sure that faces create associated edges
	if(!grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
	{
		LOG("  INFO: auto-enabling FACEOPT_AUTOGENERATE_EDGES in AdjustEdgeLength.\n");
		grid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}

	LOG("adjusting edge length\n");

//	position attachment
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);
	Grid::VertexAttachmentAccessor<APosition> aaPosRef(*pRefGrid, aPos);

//	we need an integer attachment (a helper for ObtainSimpleGrid)
	AInt aInt;
	grid.attach_to<VertexBase>(aInt);
	Grid::AttachmentAccessor<VertexBase, AInt> aaInt(grid, aInt);

//	TODO: normals shouldn't be created here but instead be passed to the method.
//	attach the vertex normals.
	ANormal aNorm;
	grid.attach_to<VertexBase>(aNorm);
	Grid::AttachmentAccessor<VertexBase, ANormal> aaNorm(grid, aNorm);
	CalculateVertexNormals(grid, aPos, aNorm);

//	assign vertex marks
	AssignFixedVertices(grid, shMarks);
	AssignCreaseVertices(grid, shMarks);

//	we need an selector that holds all edges that are candidates for a collapse
	EdgeSelector esel(grid);
	esel.enable_selection_inheritance(false);

//	sort the triangles of pRefGrid into a octree to speed-up projection performance
	PROFILE_BEGIN(octree_construction);
	SPOctree octree = CreateOctree(*pRefGrid, pRefGrid->begin<Triangle>(),
									pRefGrid->end<Triangle>(),
									15, 30, false, aPos);
	PROFILE_END();
	
	node_tree::Traverser_ProjectPoint pojectionTraverser;
	
	if(!octree.is_valid()){
		UG_LOG("  Octree creation failed in AdjustEdgeLength. Aborting.\n");
		return false;
	}
	
//	start the main iteration
	for(int iteration = 0; iteration < numIterations; ++iteration)
	{
	//	perform splits
		esel.select(grid.begin<EdgeBase>(), grid.end<EdgeBase>());
		if(!PerformSplits(grid, shMarks, esel, maxEdgeLen, aaPos, aaNorm))
			return false;

	//	perform collapses
		esel.select(grid.begin<EdgeBase>(), grid.end<EdgeBase>());
		if(!PerformCollapses(grid, shMarks, esel, minEdgeLen, aaPos, aaNorm, aaInt))
			return false;

	//	perform swaps
		esel.select(grid.begin<EdgeBase>(), grid.end<EdgeBase>());
		if(!PerformSwaps(grid, shMarks, esel, aaPos, aaNorm, aaInt))
			return false;

	//	relocate points
		LOG("  smoothing points...");
		//if(0)
		{
			PROFILE_BEGIN(smoothing_points);
			vector<vector3> vNodes;
			vector<VertexBase*> vNeighbours;
			for(int i = 0; i < 10; ++i){
				for(VertexBaseIterator iter = grid.begin<VertexBase>();
					iter != grid.end<VertexBase>(); ++iter)
				{
					VertexBase* vrt = *iter;
				//	if the vertex has marks then leave it where it is
//TODO:	crease-vertices should be moved, but should only be influenced by associated crease edges
					if(shMarks.get_subset_index(vrt) != RM_NONE)
						continue;
				
				//	collect the neighbours and project them to the plane
				//	that is defined by vrt and its normal
					vector3 v = aaPos[vrt];
					vector3 n = aaNorm[vrt];

					CollectNeighbours(vNeighbours, grid, vrt);
					vNodes.resize(vNeighbours.size());
					
					for(size_t j = 0; j < vNodes.size(); ++j)
						ProjectPointToPlane(vNodes[j], aaPos[vNeighbours[j]], v, n);
					
				//	perform point relocation
					RelocatePointBySmoothing(aaPos[vrt], v, vNodes, 5, 0.1);
				}
			}
			PROFILE_END();
		}
		LOG(" done\n");

	//	project points back on the surface
		if(projectPoints)
		//if(0)
		{
			LOG("  projecting points...");
			PROFILE_BEGIN(projecting_points);
			for(VertexBaseIterator iter = grid.vertices_begin();
				iter != grid.vertices_end(); ++iter)
			{
//TODO:	project crease vertices onto creases only! Don't project fixed vertices
				vector3 vNew;
				if(pojectionTraverser.project(aaPos[*iter], octree)){
					aaPos[*iter] = pojectionTraverser.get_closest_point();
				}
				else{
					LOG("f");
				}
/*
				if(ProjectPointToSurface(vNew, aaPos[*iter], aaNorm[*iter],
										pRefGrid->begin<Triangle>(),
										pRefGrid->end<Triangle>(), aaPosRef, true))
				{
					aaPos[*iter] = vNew;
				}
				else{
					LOG("f");
				}
*/
				
			}
			PROFILE_END();
			LOG(" done\n");
		}
	}


//	detach
	grid.detach_from<VertexBase>(aInt);
	grid.detach_from<VertexBase>(aNorm);

	return true;
}




/*
////////////////////////////////////////////////////////////////////////
//	This is an alternative version for PerformSplits.
//	It uses the Refine method.
//	It is not yet optimized for maximum speed.
//	While it performs a little less splits, overall runtime of
//	AdjustEdgeLength is not better than with the original
//	PerformSplits method.
template <class TAAPosVRT, class TAANormVRT>
bool PerformSplits(Grid& grid, SubsetHandler& shMarks, EdgeSelector& esel,
					  number maxEdgeLen, TAAPosVRT& aaPos, TAANormVRT& aaNorm)
{
	AInt aInt;
	grid.attach_to_edges(aInt);

//	compare squares
	maxEdgeLen *= maxEdgeLen;

	Selector sel(grid);
	sel.enable_autoselection(true);
	sel.select(esel.begin<EdgeBase>(), esel.end<EdgeBase>());

	LOG("  performing splits\n");
	int numSplits = 0;

	while(!sel.empty()){
	//	deselect all vertices and faces
		sel.clear_selection<VertexBase>();
		sel.clear_selection<Face>();
		sel.clear_selection<Volume>();

	//	deselect all edges that shall not be splitted
		EdgeBaseIterator iter = sel.begin<EdgeBase>();
		while(iter != sel.end<EdgeBase>()){
			EdgeBase* e = *iter;
			++iter;

		//	the higher the curvature the smaller the maxEdgeLen.
		//	minimal lenFac is 0.1
			number lenFac = CalculateLengthFac(grid, shMarks, e, aaPos);

		//	fixed edges will not be refined
			if(shMarks.get_subset_index(e) == RM_FIXED)
				sel.deselect(e);
			else if(VecDistanceSq(aaPos[e->vertex(0)], aaPos[e->vertex(1)]) < lenFac * maxEdgeLen)
				sel.deselect(e);
		}

	//	refine the grid
		Refine(grid, sel, aInt);

	//	new vertices are selected
		numSplits += sel.num<VertexBase>();

	//	re-triangulate
		Triangulate(grid, &aaPos);

	//	calculate normal for new vertices
//TODO:	be careful with crease edges
		for(VertexBaseIterator iter = sel.begin<VertexBase>();
			iter != sel.end<VertexBase>(); ++iter)
		{
			CalculateVertexNormal(aaNorm[*iter], grid, *iter, aaPos);
		}

		sel.clear_selection<VertexBase>();
		sel.clear_selection<Face>();
		sel.clear_selection<Volume>();
	}

	grid.detach_from_edges(aInt);
	LOG("  splits performed: " << numSplits << endl);
	return true;
}
*/

}//	end of namespace
