//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m01 d22

#include <fstream>
#include <queue>
#include "lib_grid/lg_base.h"
#include "common/profiler/profiler.h"
#include "simple_grid.h"
#include "edge_length_adjustment.h"
#include "lib_grid/algorithms/refinement/regular_refinement.h"
#include "common/node_tree/node_tree.h"
#include "lib_grid/algorithms/trees/octree.h"
#include "lib_grid/algorithms/callback_util.h"

using namespace std;

namespace ug
{

///	only for debugging purposes!!!
/**	Output value pairs to gnuplot...
 * \{ */
//#define SMOOTH_REMESHING__GPLOT_ENABLED
#ifdef SMOOTH_REMESHING__GPLOT_ENABLED
	typedef vector<pair<number, number> > GnuplotData;
	static GnuplotData gplotLengthFac;
	static GnuplotData gplotMinCurvature;
	static GnuplotData gplotAverageCurvature;

	void WriteGnuplotData(const char* filename, const GnuplotData& data)
	{
		ofstream out(filename);
		if(!out)
			return;

		for(size_t i = 0; i < data.size(); ++i)
			out << data[i].first << " " << data[i].second << endl;

		out.close();
	}

	#define GPLOTPOINT(dataName, x, y) dataName.push_back(make_pair<number, number>((x), (y)));
	#define GPLOTSAVE()	{WriteGnuplotData("length_fac.gplot", gplotLengthFac);\
						WriteGnuplotData("min_curvature.gplot", gplotMinCurvature);\
						WriteGnuplotData("average_curvature.gplot", gplotAverageCurvature);}
#else
//	do nothing if SMOOTH_REMESHING__GPLOT_ENABLED is false
	#define GPLOTPOINT(dataName, x, y)
	#define GPLOTSAVE()
#endif
/** \} */


/*
vector3 PNTrianglePos(const vector3& p0, const vector3& p1, const vector3& p2,
					  const vector3& n0, const vector3& n1, const vector3& n2);

vector3 PNTriangleNorm(const vector3& p0, const vector3& p1, const vector3& p2,
					   const vector3& n0, const vector3& n1, const vector3& n2);

vector3 PNCTrianglePos(const vector3& p0, const vector3& p1, const vector3& p2,
						const vector3& n0, const vector3& n1, const vector3& n2,
						const vector3& cn0, const vector3& cn1, const vector3& cn2);

vector3 PNCTriangleNorm(const vector3& p0, const vector3& p1, const vector3& p2,
						const vector3& n0, const vector3& n1, const vector3& n2,
						const vector3& cn0, const vector3& cn1, const vector3& cn2);
*/


class ILocalRemesher{
	public:
		virtual ~ILocalRemesher()	{}
		virtual void smooth_vertex(VertexBase* vrt) = 0;
		virtual VertexBase* collapse_edge(EdgeBase* edge) = 0;
		virtual VertexBase* split_edge(EdgeBase* edge) = 0;
};


class PatchRemesher : public ILocalRemesher{
	public:
		virtual ~PatchRemesher();

	///	set the grid which will be remeshed
		void set_grid(Grid& grid, APosition aPos);
		void set_crease_callbacks(Grid::vertex_traits::callback vrtCreaseCallback,
								  Grid::edge_traits::callback edgeCreaseCallback);
		void set_fixed_callbacks(Grid::vertex_traits::callback vrtFixedCallback,
								 Grid::edge_traits::callback edgeFixedCallback);

	///	Adds a new patch consisting of the given elements and associated vertices
	/**	Make sure to specify the source grid before calling this method.*/
		template <class TElemIterator>
		void add_surface_patch(TElemIterator begin, TElemIterator end);

		virtual void smooth_vertex(VertexBase* vrt);
		virtual VertexBase* collapse_edge(EdgeBase* edge);
		virtual VertexBase* split_edge(EdgeBase* edge);
		virtual EdgeBase* swap_edge(EdgeBase* edge);

		virtual vector3 vertex_position(VertexBase* vrt);
		virtual vector3 vertex_normal(VertexBase* vrt);

	protected:
		virtual void relocate_vertex(VertexBase* vrt);

	private:
		class ProjectedPoint{
			int 				patchID;
			GeometricObject* 	elem;
			vector2				barycentricCoords;
		};

		Grid	m_refGrid;
		Grid*	m_remeshGrid;
};



////////////////////////////////////////////////////////////////////////
static void AssignFixedVertices(Grid& grid, SubsetHandler& shMarks)
{	
	grid.begin_marking();
	
//	mark all vertices contained in a crease edge,
//	and which are not regular crease-vertices as fixed
	for(EdgeBaseIterator iter = shMarks.begin<EdgeBase>(REM_CREASE);
		iter != shMarks.end<EdgeBase>(REM_CREASE); ++iter)
	{
		EdgeBase* e = *iter;
		for(size_t i = 0; i < 2; ++i){
			VertexBase* vrt = e->vertex(i);
			if(!grid.is_marked(vrt)){
				grid.mark(vrt);
				int counter = 0;
				for(Grid::AssociatedEdgeIterator nbIter = grid.associated_edges_begin(vrt);
					nbIter != grid.associated_edges_end(vrt); ++nbIter)
				{
					if(shMarks.get_subset_index(*nbIter) != REM_NONE)
						++counter;
				}
				
				if(counter != 2)
					shMarks.assign_subset(vrt, REM_FIXED);
			}
		}
	}
	
	grid.end_marking();
	
//	mark all vertices that lie on a fixed edge as fixed
	for(EdgeBaseIterator iter = shMarks.begin<EdgeBase>(REM_FIXED);
		iter != shMarks.end<EdgeBase>(REM_FIXED); ++iter)
	{
		EdgeBase* e = *iter;
		shMarks.assign_subset(e->vertex(0), REM_FIXED);
		shMarks.assign_subset(e->vertex(1), REM_FIXED);
	}
}

static void AssignCreaseVertices(Grid& grid, SubsetHandler& shMarks)
{
//	mark all vertices that lie on a crease and which are not fixed
//	as crease vertices.
	if(shMarks.num_subsets() <= REM_CREASE)
		return;

	for(EdgeBaseIterator iter = shMarks.begin<EdgeBase>(REM_CREASE);
		iter != shMarks.end<EdgeBase>(REM_CREASE); ++iter)
	{
		EdgeBase* e = *iter;
		for(uint i = 0; i < 2; ++i)
			if(shMarks.get_subset_index(e->vertex(i)) != REM_FIXED)
				shMarks.assign_subset(e->vertex(i), REM_CREASE);
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
	Grid::AssociatedFaceIterator iterEnd = grid.associated_faces_end(vrt);
	for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(vrt);
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

//todo:	Think about converting to radiants.
	//minDot = acos(minDot) / (0.5 * PI);
	//...
//	done
	GPLOTPOINT(gplotMinCurvature, 0, minDot);
	return minDot;	
}							

////////////////////////////////////////////////////////////////////////
template <class TAAPosVRT>
number CalculateAverageCurvature(Grid& grid, SubsetHandler& shMarks,
								EdgeBase* e, TAAPosVRT& aaPos)
{
	number avCurv = 0.5 * (CalculateMinCurvature(grid, shMarks, e->vertex(0), aaPos)
				+ CalculateMinCurvature(grid, shMarks, e->vertex(1), aaPos));
	GPLOTPOINT(gplotAverageCurvature, 0.5, avCurv);
	return avCurv;
}

////////////////////////////////////////////////////////////////////////
template <class TAAPosVRT>
number CalculateLengthFac(Grid& grid, SubsetHandler& shMarks,
								EdgeBase* e, TAAPosVRT& aaPos)
{
	number lenFac = CalculateAverageCurvature(grid, shMarks, e, aaPos);
	lenFac = (lenFac - 0.95) / 0.05;
	lenFac = max(number(0.25), lenFac);
	GPLOTPOINT(gplotLengthFac, 0.5, lenFac);
	return lenFac;
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
		if(pshMarks->get_subset_index(e) == REM_FIXED ||
			pshMarks->get_subset_index(e) == REM_CREASE)
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
	//number approxDeg = GeometricApproximationDegree(sg);
	number shapeDeg = ShapeQualityDegree(sg);

	number smoothDeg = VecDot(sg.triangleNormals[0], sg.triangleNormals[1]);

//	this is a new test. the idea is that each edge should be orthogonal to
//	the normals of its endpoints - at least in a perfectly smooth surface.

	number approxDeg;
	number newApproxDeg;
	{
		vector3 v;
		VecSubtract(v, sg.vertices[1], sg.vertices[0]);
		VecNormalize(v, v);
		approxDeg = 1. - 0.5*(fabs(VecDot(v, sg.vertexNormals[0])) +
							fabs(VecDot(v, sg.vertexNormals[1])));
						
		VecSubtract(v, sg.vertices[3], sg.vertices[2]);
		VecNormalize(v, v);
		newApproxDeg = 1. - 0.5*(fabs(VecDot(v, sg.vertexNormals[2])) +
								fabs(VecDot(v, sg.vertexNormals[3])));
	}

//	perform a swap on the simple grid
	if(!SwapEdge(sg))
	{
		LOG("swap edge failed...\n");
		return false;
	}

//	calculate new geometric-approximation-degree and triangle quality
	//number newApproxDeg = GeometricApproximationDegree(sg);
	number newShapeDeg = ShapeQualityDegree(sg);
	number newSmoothDeg = VecDot(sg.triangleNormals[0], sg.triangleNormals[1]);

//	neither the shapeDeg nor the approxDeg may get too bad.
	if((newApproxDeg < 0.5 * approxDeg) || (newShapeDeg < 0.5 * shapeDeg))
		return false;
	
//	make sure that the swap does not destroy the smoothness
	if(newSmoothDeg < 0.1 * smoothDeg)
		return false;

	if(0.2 * (newApproxDeg - approxDeg) + 0.8 * (newShapeDeg - shapeDeg) > 0)
	//if(newShapeDeg > shapeDeg)
	//if(newApproxDeg > approxDeg)
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
	int maxNumSwaps = esel.num<EdgeBase>() * 2;
		
	while(!esel.empty())
	{
		EdgeBase* e = *esel.begin<EdgeBase>();
		esel.deselect(e);

		if(TrySwap(grid, e, aaPos, aaNorm, aaInt, &shMarks, &esel)){
			++numSwaps;
			if(numSwaps > maxNumSwaps){
				UG_LOG("  aborting since maxNumSwaps was reached...");
				esel.clear();
			}
		}
	}
	LOG("  swaps performed: " << numSwaps << endl);
	
	return true;
}

/**	returns the resulting vertex or NULL, if no collapse was performed.*/
template <class TAAPosVRT, class TAANormVRT, class TAAIntVRT>
VertexBase* TryCollapse(Grid& grid, EdgeBase* e,
				TAAPosVRT& aaPos, TAANormVRT& aaNorm, 
				TAAIntVRT& aaInt, SubsetHandler* pshMarks = NULL,
				EdgeSelector* pCandidates = NULL)
{
	if(pshMarks)
	{
		SubsetHandler& shMarks = *pshMarks;
	//	collapses are not allowed for fixed edges
		if(shMarks.get_subset_index(e) == REM_FIXED)
			return NULL;
			
	//	if both endpoints of are fixed vertices then
	//	we may not collapse
		int vrtSI[2];
		vrtSI[0] = shMarks.get_subset_index(e->vertex(0));
		vrtSI[1] = shMarks.get_subset_index(e->vertex(1));

		if((vrtSI[0] == REM_FIXED) && (vrtSI[1] == REM_FIXED))
			return NULL;

	//	if both endpoints are somehow marked, e has to be a
	//	crease edge
		if((vrtSI[0] != REM_NONE) && (vrtSI[1] != REM_NONE)
			&&	(shMarks.get_subset_index(e) != REM_CREASE))
			return NULL;
	}

//	check whether the edge can be collapsed
	if(EdgeCollapseIsValid(grid, e))
	{
	//	test the collapse using a simple-grid
		SimpleGrid sgSrc;
		if(!ObtainSimpleGrid(sgSrc, grid, e->vertex(0), e->vertex(1), 1,
							aaPos, aaNorm, aaInt))
		{
			LOG("ObtainSimpleGrid failed. ignoring edge...\n");
			return NULL;
		}
		
	//	calculate geometric-approximation-degree and triangle quality
		number approxDeg = GeometricApproximationDegree(sgSrc);
		number shapeDeg = ShapeQualityDegree(sgSrc);

	//	perform a collapse on the simple grid
		SimpleGrid sgDest;
		if(!ObtainSimpleGrid_CollapseEdge(sgDest, grid, e,
									1, aaPos, aaNorm, aaInt))
		{
			LOG("collapse edge failed...\n");
			return NULL;
		}

	//	get the positions of the old endpoints
		static const int numTestPositions = 3;
		int newInd = 0;
		vector3 v[numTestPositions];
		v[0] = aaPos[e->vertex(0)];
		v[1] = aaPos[e->vertex(1)];
		v[2] = sgDest.vertices[newInd];
		
	//	we'll compare 3 approximation degrees and three shape degrees
		number newApproxDeg[numTestPositions];
		number newShapeDeg[numTestPositions];
		
	//	check which position is the best
		int bestIndex = -1;
	//	the vertex subset index is used to support marks (crease and fixed vertices)
		int vrtSI[2];
		vrtSI[0] = vrtSI[1] = REM_NONE;

		if(pshMarks)
		{
			vrtSI[0] = pshMarks->get_subset_index(e->vertex(0));
			vrtSI[1] = pshMarks->get_subset_index(e->vertex(1));

			if((vrtSI[0] == REM_FIXED) || ((vrtSI[0] != REM_NONE) && (vrtSI[1] == REM_NONE))){
				bestIndex = 0;
				sgDest.vertices[newInd] = v[0];
				newApproxDeg[0] = GeometricApproximationDegree(sgDest);
				newShapeDeg[0] = ShapeQualityDegree(sgDest);
			}
			else if((vrtSI[1] == REM_FIXED) || ((vrtSI[1] != REM_NONE) && (vrtSI[0] == REM_NONE))){
				bestIndex = 1;
				sgDest.vertices[newInd] = v[1];
				newApproxDeg[1] = GeometricApproximationDegree(sgDest);
				newShapeDeg[1] = ShapeQualityDegree(sgDest);
			}
		}

		if(bestIndex == -1){
		//	check all three approximation degrees
			for(int i = 0; i < numTestPositions; ++i){
			//	we'll compute all qualities with the averaged normal
				sgDest.vertices[newInd] = v[i];
				CalculateTriangleNormals(sgDest);
				newApproxDeg[i] = GeometricApproximationDegree(sgDest);
				newShapeDeg[i] = ShapeQualityDegree(sgDest);
			}
		//	get the best one
			bestIndex = 0;
		/*
			for(int i = 1; i < numTestPositions; ++i){
				if(newApproxDeg[i] > newApproxDeg[bestIndex])
					bestIndex = i;
			}
		*/
			for(int i = 1; i < numTestPositions; ++i){
				if(newShapeDeg[i] > newShapeDeg[bestIndex])
					bestIndex = i;
			}
		}

	//	if the shape-degree of the collapsed region is too bad, we'll skip the collapse
		if(newShapeDeg[bestIndex] < 0.5 * shapeDeg)
			return NULL;
/*
	//	the approximation degree is only interesting if both endpoints of the
	//	edge are regular surface vertices
		bool regularNeighbourhood = IsRegularSurfaceVertex(grid, e->vertex(0)) &&	
									IsRegularSurfaceVertex(grid, e->vertex(1));
*/

	//	if the best approximation degree is not too bad, we'll perform the collapse
		if(/*!regularNeighbourhood || */(newApproxDeg[bestIndex] > 0.95 * approxDeg))
		{						
		//	pick one of the endpoints to be the one that resides
		//	This has to be done with care, since the residing vertex
		//	determines which edges will be removed and which will remain
		//	(this is important, since we want crease edges to remain!).
		//	If only one of the endpoints lies on a crease or is fixed, the
		//	decision is easy (the marked one has to remain).
		//	If both are marked, we have to investigate the triangles which
		//	are adjacent to the collapsed edge (quadrilaterals don't make
		//	problems here).
		//	- If all three corner vertices are marked we have to distinguish:
		//		(a) all three edges are marked -> abort
		//		(b) two edges are marked -> the vertex connecting both remains.
		//		(c) one edge is marked -> the chosen vertex doesn't affect edge-marks.
		//	- If only the two are marked they don't affect the edge-marks.
		//	
		//	Even more care has to be taken if a fixed vertex is involved.
		//	If case (b) applies and the creas-vertex is has to remain, it
		//	has to be moved to the position of the fixed-vertex and has to
		//	be marked fixed itself.
					
		//	choose the vertex that shall remain.
			VertexBase* vrt = e->vertex(0);

			if(vrtSI[0] != REM_FIXED && vrtSI[1] != REM_NONE)
			{
				vrt = e->vertex(1);
			}

			if(vrtSI[0] != REM_NONE && vrtSI[1] != REM_NONE){
			//	both are marked. Things are getting complicated now!
			//	get adjacent faces of e
				vector<Face*> vFaces;
				vFaces.reserve(2);
				CollectFaces(vFaces, grid, e);
			
				vector<EdgeBase*> vEdges;
				vEdges.reserve(4);
				for(size_t i = 0; i < vFaces.size(); ++i){
					Face* f = vFaces[i];
				//	only triangles are interesting
					if(f->num_edges() == 3){
					//	get the number of marked edges
						CollectEdges(vEdges, grid, f);
					//	count the marked edges
						int numMarked = 0;
						for(size_t j = 0; j < vEdges.size(); ++j){
						//	note that pshMarks exists since vrtSI != REM_NONE
							if(pshMarks->get_subset_index(vEdges[j]) != REM_NONE)
								++numMarked;
						}
					
						if(numMarked == 3){
							return NULL;	//	case (a) applies
						}
						else if(numMarked == 2){
						//	check which of the vrts is connected to two
						//	marked edges of the triangle
							for(size_t j = 0; j < 2; ++j){
								int numMarked = 0;
								for(size_t k = 0; k < vEdges.size(); ++k){
									if(pshMarks->get_subset_index(vEdges[k]) != REM_NONE){
										if(EdgeContains(vEdges[k], e->vertex(j)))
											++numMarked;
									}
								}
							//	if numMarked == 2 we found it
								if(numMarked == 2){
								//	the connected edge has to be marked as a crease
									EdgeBase* ce = GetConnectedEdge(grid, e->vertex(j), f);
									if(ce)
										pshMarks->assign_subset(ce, REM_CREASE);
								//	we're done. break
									break;
								}
							}
						}
						else{
						}
					}
				}
			}

		//	collapse the edge
			CollapseEdge(grid, e, vrt);
			
		//	assign best position
			aaPos[vrt] = v[bestIndex];
		//	assign the normal
			aaNorm[vrt] = sgDest.vertexNormals[newInd];
/*
			if(pCandidates){
//TODO: all edges that belong to associated faces are new candidates.						
			//	associated edges of vrt are possible new candidates
				pCandidates->select(grid.associated_edges_begin(vrt),
							grid.associated_edges_end(vrt));
			}*/

			return vrt;
		}
	}

	return NULL;
}

////////////////////////////////////////////////////////////////////////
template <class TAAPosVRT, class TAANormVRT, class TAAIntVRT>
bool PerformCollapses(Grid& grid, SubsetHandler& shMarks, EdgeSelector& esel,
					  number minEdgeLen, TAAPosVRT& aaPos, TAANormVRT& aaNorm,
					  TAAIntVRT& aaInt, bool adaptive = true)
{	
	PROFILE_FUNC();
	LOG("  performing collapses\n");
	vector<EdgeBase*>	edges;
	int numCollapses = 0;
//	compare squares
	minEdgeLen *= minEdgeLen;

	while(!esel.empty())
	{
		EdgeBase* e = *esel.begin<EdgeBase>();
		esel.deselect(e);
		
	//	the higher the curvature the smaller the maxEdgeLen.
	//	minimal lenFac is 0.1
		number lenFac = 1.0;
		if(adaptive)
			lenFac = CalculateLengthFac(grid, shMarks, e, aaPos);

	//	check whether the edge is short enough
		if(VecDistanceSq(aaPos[e->vertex(0)], aaPos[e->vertex(1)]) < lenFac * minEdgeLen)
		{
			VertexBase* vrt = TryCollapse(grid, e, aaPos, aaNorm, aaInt, &shMarks, &esel);
			if(vrt){
				++numCollapses;

			//	we'll deselect associated edges of vrt, to avoid cascade-collapses
				CollectAssociated(edges, grid, vrt);
				esel.deselect(edges.begin(), edges.end());
			}
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
		if(pshMarks->get_subset_index(e) == REM_FIXED)
			return false;
		else if(pshMarks->get_subset_index(e) == REM_CREASE)
			bCreaseEdge = true;
	}

//	get the center of the edges
	vector3 vCenter = CalculateCenter(e, aaPos);

//	the new normal
	vector3 n;
	if(bCreaseEdge){		
	//	interpolating the normal can cause severe problems at highly
	//	irregular vertices or if one vertecs lies on a very
	//	sharp edge (the normals of the endpoints thus point
	//	in different directions.)
//		VecAdd(n, aaNorm[e->vertex(0)], aaNorm[e->vertex(1)]);
//		VecNormalize(n, n);
		CalculateNormal(n, grid, e, aaPos);
	}

//	split the edge
	Vertex* vrt = SplitEdge<Vertex>(grid, e, false);

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
					  number maxEdgeLen, TAAPosVRT& aaPos, TAANormVRT& aaNorm,
					  bool adaptive = true)
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
		number lenFac = 1.0;
		if(adaptive)
			lenFac = CalculateLengthFac(grid, shMarks, e, aaPos);

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
	if(vNodes.empty())
		return;

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
//	FixBadTriangles
template <class TAAPosVRT, class TAANormVRT>
bool FixBadTriangles(Grid& grid, SubsetHandler& shMarks, EdgeSelector& esel,
					TAAPosVRT& aaPos, TAANormVRT& aaNorm,
					number qualityThreshold)
{
	LOG("  fixing bad triangles... ");
//	iterate through marked edges and check whether adjacent triangles are
//	badly shaped. If this is the case, then split the non-marked
//	edges of the triangle.
//	store newly inserted vertices and smooth them at the end of the algo
	vector<VertexBase*> vNewVrts;
	vector<Face*> vFaces;
	vector<EdgeBase*> vEdges;
	
//	we wont assign this selector to the grid until it is clear that
//	we'll need it.
	Selector sel;
	
	
	EdgeBaseIterator iter = esel.begin<EdgeBase>();
	while(iter != esel.end<EdgeBase>())
	{
	//	store edge and increase iterator immediatly
		EdgeBase* e = *iter;
		++iter;
		
	//	get the adjacent faces
		CollectFaces(vFaces, grid, e);
		
	//	check whether one of them is degenerated.
		for(size_t i = 0; i < vFaces.size(); ++i){
			Face* f = vFaces[i];
			number q = FaceQuality(f, aaPos);
		//	if the quality is smaller than threshold, we have to
		//	do something
			if(q < qualityThreshold){
			//	make sure that the selector is connected to the grid
				if(sel.grid() == NULL)
					sel.assign_grid(grid);
					
			//	get associated edges and mark them for refinement
				CollectEdges(vEdges, grid, f);

				for(size_t j = 0; j < vEdges.size(); ++j){
					if(shMarks.get_subset_index(vEdges[j]) == REM_NONE)
						sel.select(vEdges[j]);
				}
			}
		}
	}
	
//	if edges have been selected, we'll now call refine.
	if(sel.grid() != NULL){
		if(sel.num<EdgeBase>() > 0){
			if(Refine(grid, sel)){
				LOG(sel.num<VertexBase>() << " new vertices... ");
			//	retriangulate surface
				if(grid.num<Quadrilateral>() > 0)
					Triangulate(grid, grid.begin<Quadrilateral>(), grid.end<Quadrilateral>());
				
			//	calculate normals, then
			//	smooth new vertices (all vertices selected in sel).
				vector<VertexBase*> vNeighbours;
				vector<vector3> vNodes;
				
			//	calculate normals
				for(VertexBaseIterator iter = sel.begin<VertexBase>();
						iter != sel.end<VertexBase>(); ++iter)
				{
				//	collect neighbour nodes
					VertexBase* vrt = *iter;
					CollectNeighbors(vNeighbours, grid, vrt);
					
				//	sum their normals and interpolate it
				//	make sure to only add valid normals
					vector3 n(0, 0, 0);
					for(size_t i = 0; i < vNeighbours.size(); ++i){
						if(!sel.is_selected(vNeighbours[i]))
							VecAdd(n, n, aaNorm[vNeighbours[i]]);
					}
					VecNormalize(aaNorm[vrt], n);
				}
				
			//	repeat smoothing.
				for(size_t i = 0; i < 5; ++i){
					for(VertexBaseIterator iter = sel.begin<VertexBase>();
						iter != sel.end<VertexBase>(); ++iter)
					{
						VertexBase* vrt = *iter;

					//	collect the neighbours and project them to the plane
					//	that is defined by vrt and its normal
						vector3 v = aaPos[vrt];
						vector3 n = aaNorm[vrt];

						CollectNeighbors(vNeighbours, grid, vrt);
						vNodes.resize(vNeighbours.size());
						
						for(size_t j = 0; j < vNodes.size(); ++j)
							ProjectPointToPlane(vNodes[j], aaPos[vNeighbours[j]], v, n);
						
					//	perform point relocation
						RelocatePointBySmoothing(aaPos[vrt], v, vNodes, 5, 0.05);
					}
				}
			}
			else{
				LOG("refine failed!\n");
				return false;
			}
		}
	}
	
	LOG("done\n");
	return true;
}

////////////////////////////////////////////////////////////////////////
//	PerformSmoothing
template <class TAAPosVRT, class TAANormVRT>
void PerformSmoothing(Grid& grid, SubsetHandler& shMarks,
					TAAPosVRT& aaPos, TAANormVRT& aaNorm,
					size_t numIterations, number stepSize)
{
	vector<vector3> vNodes;
	vector<VertexBase*> vNeighbours;
	for(size_t i = 0; i < numIterations; ++i){
		CalculateVertexNormals(grid, aaPos, aaNorm);
		for(VertexBaseIterator iter = grid.begin<VertexBase>();
			iter != grid.end<VertexBase>(); ++iter)
		{
			VertexBase* vrt = *iter;
		//	if the vertex is fixed then leave it where it is.
			if(shMarks.get_subset_index(vrt) == REM_FIXED)
				continue;
/*
if(shMarks.get_subset_index(vrt) == REM_CREASE)
	continue;
*/
		//	collect the neighbours and project them to the plane
		//	that is defined by vrt and its normal
			vector3 v = aaPos[vrt];
			vector3 n = aaNorm[vrt];

			bool bProjectPointsToPlane = true;
			
			if(shMarks.get_subset_index(vrt) == REM_CREASE){
				CollectNeighbors(vNeighbours, grid, vrt, NHT_EDGE_NEIGHBORS,
									IsInSubset(shMarks, REM_CREASE));

			//	we have to choose a special normal
				if(vNeighbours.size() != 2){
					UG_LOG("n"<<vNeighbours.size());
					continue;
				}
				else{
					vector3 v1, v2;
					VecSubtract(v1, v, aaPos[vNeighbours[0]]);
					VecSubtract(v2, v, aaPos[vNeighbours[1]]);
					VecNormalize(v1, v1);
					VecNormalize(v2, v2);
					VecAdd(n, v1, v2);
					if(VecLengthSq(n) > SMALL)
						VecNormalize(n, n);
					else {
					//	both edges have the same direction.
					//	don't project normals
						bProjectPointsToPlane = false;
					}
				}
			}
			else{
				CollectNeighbors(vNeighbours, grid, vrt);
			}
			
			vNodes.resize(vNeighbours.size());

			if(bProjectPointsToPlane){
				for(size_t j = 0; j < vNodes.size(); ++j)
					ProjectPointToPlane(vNodes[j], aaPos[vNeighbours[j]], v, n);
			}
			else{
				for(size_t j = 0; j < vNodes.size(); ++j)
					vNodes[j] = aaPos[vNeighbours[j]];
			}
		//	perform point relocation
			RelocatePointBySmoothing(aaPos[vrt], v, vNodes, 5, stepSize);
		}
	}
}

/**	Make sure that elements in gridOut directly correspond to
 *	elements in gridIn*/
template <class TGeomObj>
void CopySelectionStatus(Selector& selOut, Grid& gridOut,
						 Selector& selIn, Grid& gridIn)
{
	typedef typename geometry_traits<TGeomObj>::iterator GeomObjIter;
	GeomObjIter endOut = gridOut.end<TGeomObj>();
	GeomObjIter endIn = gridIn.end<TGeomObj>();
	GeomObjIter iterOut = gridOut.begin<TGeomObj>();
	GeomObjIter iterIn = gridIn.begin<TGeomObj>();
	
	for(; (iterOut != endOut) && (iterIn != endIn); ++iterOut, ++iterIn)
	{
		if(selIn.is_selected(*iterIn))
			selOut.select(*iterOut);
	}
}

////////////////////////////////////////////////////////////////////////
bool AdjustEdgeLength(Grid& grid, SubsetHandler& shMarks,
					  number minEdgeLen, number maxEdgeLen, int numIterations,
					  bool projectPoints, bool adaptive)
{
	PROFILE_FUNC();

//TODO:	replace this by a template parameter
	APosition aPos = aPosition;

//	we have to make sure that the mesh consist of triangles only,
//	since the algorithm would produce bad results if not.
	if(grid.num<Quadrilateral>() > 0)
	{
		UG_LOG("  INFO: grid contains quadrilaterals. Converting to triangles...\n");
				
//TODO:	not gridIn but a copy of gridIn (pRefGrid) should be triangulated.
		Triangulate(grid, grid.begin<Quadrilateral>(),
					grid.end<Quadrilateral>());
	}
		
//	make sure that grid and pRefGrid have position-attachments
	if(!grid.has_vertex_attachment(aPos)){
		UG_LOG("  vertex-position-attachment missing in AdjustEdgeLength. Aborting.\n");
		return false;
	}
	
//	make sure that faces create associated edges
	if(!grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
	{
		LOG("  INFO: auto-enabling FACEOPT_AUTOGENERATE_EDGES in AdjustEdgeLength.\n");
		grid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}

//	position attachment
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);

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

//	we need an selector that holds all edges that are candidates for a operation
	EdgeSelector esel(grid);
	esel.enable_selection_inheritance(false);
	
//	sort the triangles of pRefGrid into a octree to speed-up projection performance
	SPOctree octree;
	node_tree::Traverser_ProjectPoint pojectionTraverser;
	if(projectPoints){
		//PROFILE_BEGIN(octree_construction);
		octree = CreateOctree(grid, grid.begin<Triangle>(),
									grid.end<Triangle>(),
									10, 30, false, aPos);

		//PROFILE_END();
		if(!octree.valid()){
			UG_LOG("  Octree creation failed in AdjustEdgeLength. Aborting.\n");
			return false;
		}
	}
	
//	start the main iteration
	for(int iteration = 0; iteration < numIterations; ++iteration)
	{
	//	perform splits
		esel.select(grid.begin<EdgeBase>(), grid.end<EdgeBase>());
		if(!PerformSplits(grid, shMarks, esel, maxEdgeLen, aaPos, aaNorm, adaptive))
			return false;

	//	perform collapses
		esel.select(grid.begin<EdgeBase>(), grid.end<EdgeBase>());
		if(!PerformCollapses(grid, shMarks, esel, minEdgeLen, aaPos, aaNorm, aaInt, adaptive))
			return false;

	//	perform swaps
		esel.select(grid.begin<EdgeBase>(), grid.end<EdgeBase>());
		if(!PerformSwaps(grid, shMarks, esel, aaPos, aaNorm, aaInt))
			return false;

/*
//	This is commented out, since it didn't help with the problems encountered
//	in the geometries at that time.
//	The algorithm however works and may prove useful in the future.
	//	fix bad triangles
	//	adjacent to crease-edges badly shaped triangles may occur,
	//	which have to be treated somehow.
		esel.clear();
		esel.select(shMarks.begin<EdgeBase>(REM_CREASE), shMarks.end<EdgeBase>(REM_CREASE));
		esel.select(shMarks.begin<EdgeBase>(REM_FIXED), shMarks.end<EdgeBase>(REM_FIXED));
		FixBadTriangles(grid, shMarks, esel, aaPos, aaNorm, 0.1);
*/
	//	relocate points
		/*LOG("  smoothing points...");
		PerformSmoothing(grid, shMarks, aaPos, aaNorm, 10, 0.1);
		LOG(" done\n");*/

		LOG("  updating normals...\n");
		CalculateVertexNormals(grid, aPos, aNorm);

	//	project points back on the surface
		if(projectPoints)
		{
			LOG("  projecting points...");
			//PROFILE_BEGIN(projecting_points);
			for(VertexBaseIterator iter = grid.vertices_begin();
				iter != grid.vertices_end(); ++iter)
			{
//TODO:	project crease vertices onto creases only! Don't project fixed vertices
				if(shMarks.get_subset_index(*iter) != REM_FIXED){
					vector3 vNew;
					if(pojectionTraverser.project(aaPos[*iter], octree/*, &aaNorm[*iter]*/)){
						aaPos[*iter] = pojectionTraverser.get_closest_point();
					}
					else{
						LOG("f");
					}
				}
			}
			//PROFILE_END();
			LOG(" done\n");
		}

		if(iteration < numIterations - 1){
			LOG("  updating normals...");
			CalculateVertexNormals(grid, aPos, aNorm);
		}
	}


//	detach
	grid.detach_from<VertexBase>(aInt);
	grid.detach_from<VertexBase>(aNorm);

	GPLOTSAVE();
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
			if(shMarks.get_subset_index(e) == REM_FIXED)
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
