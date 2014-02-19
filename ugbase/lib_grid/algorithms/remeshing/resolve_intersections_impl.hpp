// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Jun 21, 2013 (d,m,y)

#ifndef __H__UG__resolve_intersections_impl__
#define __H__UG__resolve_intersections_impl__

#include "resolve_intersections.h"

namespace ug{

template <class TAAPosVRT>
Vertex* ResolveVertexEdgeIntersection(Grid& grid, Vertex* v,
										   Edge* e, TAAPosVRT& aaPos,
										   number snapThreshold)
{
	typedef typename TAAPosVRT::ValueType vector_t;

	number snapThresholdSq = snapThreshold * snapThreshold;

//	make sure that the vertex is not an endpoint of e
	if(EdgeContains(e, v))
		return NULL;

//	we have to make sure that v and e are not connected by a face.
//	This could lead to infinite recursion
/*
	vector<Face*> faces;
	CollectFaces(faces, grid, e);
	for(size_t i = 0; i < faces.size(); ++i){
		if(FaceContains(faces[i], v))
			return NULL;
	}
*/
//	project the vertex to the line defined by the edge
	vector_t p;
	number t = DropAPerpendicular(p, aaPos[v], aaPos[e->vertex(0)],
								  aaPos[e->vertex(1)]);

	if((t >= 0) && (t <= 1.)){
		if(VecDistanceSq(p, aaPos[v]) < snapThresholdSq){
		//	to make sure that no double edges may occur, we'll use MergeVertices
			RegularVertex* nVrt = SplitEdge<RegularVertex>(grid, grid, e);
			aaPos[v] = p;
			MergeVertices(grid, v, nVrt);
			return v;
/*
		//	insert the vertex into the edge
			CreateEdgeSplitGeometry(grid, grid, e, v);
			grid.erase(e);
			return v;
*/
		}
	}
	return NULL;
}

/**
 * No support for volumes in the current version.
 * \todo Instead of manually refining the face, an external function SplitFace
 *		 should be used, which can take care of volumes, too.
 */
template <class TAAPosVRT>
bool ResolveVertexFaceIntersection(Grid& grid, Vertex* v,
								   Face* f, TAAPosVRT& aaPos,
								   number snapThreshold)
{
	using namespace std;
	typedef typename TAAPosVRT::ValueType vector_t;

	number snapThresholdSq = snapThreshold * snapThreshold;

//	make sure that the vertex is not a corner of f
	if(FaceContains(f, v))
		return false;

//	calculate the normal
	vector_t n;
	CalculateNormal(n, f, aaPos);

//	project the vertex to the plane defined by the face
	vector_t p;
	ProjectPointToPlane(p, aaPos[v], aaPos[f->vertex(0)], n);

//	check whether the distance is fine
	if(VecDistanceSq(p, aaPos[v]) < snapThresholdSq){
		bool refined = false;
		vector<Face*> newFaces;
		Vertex* newFaceVrt = NULL;
		Vertex* nVrt = NULL;
		vector_t pi;
	//	now we have to check whether the projection lies in the face
		if(f->num_vertices() == 3){
			if(RayTriangleIntersection(pi, aaPos[f->vertex(0)], aaPos[f->vertex(1)],
										aaPos[f->vertex(2)], p, n))
			{
			//	ok we have to insert the vertex
			//	we'll create a temporary new vertex, which will then be merged with v.
			//	This is important, since we can avoid double elements this way.
				nVrt = *grid.create<RegularVertex>();
				Vertex* newEdgeVrts[3] = {NULL, NULL, NULL};
				refined = f->refine(newFaces, &newFaceVrt, newEdgeVrts, nVrt);
			}
		}
		else if(f->num_vertices() == 4){
			bool success = false;
			if(RayTriangleIntersection(pi, aaPos[f->vertex(0)], aaPos[f->vertex(1)],
										aaPos[f->vertex(2)], p, n))
			{
				success = true;
			}
			else if(RayTriangleIntersection(pi, aaPos[f->vertex(0)], aaPos[f->vertex(2)],
											aaPos[f->vertex(3)], p, n))
			{
				success = true;
			}

			if(success){
			//	ok we have to insert the vertex
			//	we'll create a temporary new vertex, which will then be merged with v.
			//	This is important, since we can avoid double elements this way.
				nVrt = *grid.create<RegularVertex>();
				Vertex* newEdgeVrts[4] = {NULL, NULL, NULL, NULL};
				refined = f->refine(newFaces, &newFaceVrt, newEdgeVrts, nVrt);
			}
		}

		if(refined){
		//	adjust position
			aaPos[v] = pi;
		//	register the new faces and erase the old one
			for(size_t i = 0; i < newFaces.size(); ++i)
				grid.register_element(newFaces[i], f);
			grid.erase(f);

		//	to make sure that no double edges may occur, we'll use MergeVertices
			MergeVertices(grid, v, nVrt);

			return true;
		}
	}

	return false;
}

/**
 * This method does not resolve intersections between close, parallel edges or
 * between degenerate edges. You can treat such cases with
 * ReolveVertexEdgeIntersection.
 */
template <class TAAPosVRT>
Vertex* ResolveEdgeEdgeIntersection(Grid& grid, Edge* e1, Edge* e2,
										TAAPosVRT& aaPos, number snapThreshold)
{
	typedef typename TAAPosVRT::ValueType vector_t;

//	check whether one edge contains a vertex of another edge
	if(EdgeContains(e1, e2->vertex(0)) || EdgeContains(e1, e2->vertex(1)))
		return NULL;

	number snapThresholdSq = snapThreshold * snapThreshold;

	number t1, t2;
	if(LineLineProjection(t1, t2, aaPos[e1->vertex(0)], aaPos[e1->vertex(1)],
						  aaPos[e2->vertex(0)], aaPos[e2->vertex(1)]))
	{
	//	calculate the positions
		vector_t v1, v2;
		VecScaleAdd(v1, (1. - t1), aaPos[e1->vertex(0)], t1, aaPos[e1->vertex(1)]);
		VecScaleAdd(v2, (1. - t2), aaPos[e2->vertex(0)], t2, aaPos[e2->vertex(1)]);

	//	check whether the points are close to each other
		if(VecDistanceSq(v1, v2) < snapThresholdSq){
		//	calculate center
			vector_t p;
			VecScaleAdd(p, 0.5, v1, 0.5, v2);

		//	to make sure that no double edges may occur, we'll use MergeVertices
			RegularVertex* nVrt1 = SplitEdge<RegularVertex>(grid, grid, e1);
			RegularVertex* nVrt2 = SplitEdge<RegularVertex>(grid, grid, e2);
			aaPos[nVrt1] = p;
			MergeVertices(grid, nVrt1, nVrt2);

			return nVrt1;
		/*
		//	create a new vertex and split both edges using it
			RegularVertex* vrt = *grid.create<RegularVertex>();
			aaPos[vrt] = p;
			CreateEdgeSplitGeometry(grid, grid, e1, vrt);
			CreateEdgeSplitGeometry(grid, grid, e2, vrt);
			grid.erase(e1);
			grid.erase(e2);
		*/
		}
		else{
		/*
			LOG("distance check failed at: " << v1 << ", " << v2 << endl);
			UG_LOG("edges with vertices: " << aaPos[e1->vertex(0)] << aaPos[e1->vertex(1)] << endl);
			UG_LOG("                     " << aaPos[e2->vertex(0)] << aaPos[e2->vertex(1)]);
		*/
		}
	}
	return NULL;
}

/**
 * No support for volumes in the current version.
 * \todo Instead of manually refining the face, an external function SplitFace
 *		 should be used, which can take care of volume, too.
 */
template <class TAAPosVRT>
bool ResolveEdgeFaceIntersection(Grid& grid, Edge* e, Face* f,
								 TAAPosVRT& aaPos, number snapThreshold)
{
	using namespace std;
	typedef typename TAAPosVRT::ValueType vector_t;

//	check whether one edge contains a vertex of another edge
	if(FaceContains(f, e->vertex(0)) || FaceContains(f, e->vertex(1)))
		return false;

	number snapThresholdSq = snapThreshold * snapThreshold;

	vector_t dir;
	VecSubtract(dir, aaPos[e->vertex(1)], aaPos[e->vertex(0)]);

	vector_t p;
	number t1, t2, s;
	bool refined = false;
	vector<Face*> newFaces;
	Vertex* newFaceVrt = NULL;
	Vertex* vrt = NULL;
	if(f->num_vertices() == 3){
		if(RayTriangleIntersection(p, t1, t2, s, aaPos[f->vertex(0)], aaPos[f->vertex(1)],
									aaPos[f->vertex(2)], aaPos[e->vertex(0)], dir))
		{
			if((s >= 0) && (s <= 1.)){
			//	split the face
				vrt = *grid.create<RegularVertex>();
				Vertex* newEdgeVrts[3] = {NULL, NULL, NULL};
				refined = f->refine(newFaces, &newFaceVrt, newEdgeVrts, vrt);
			}
		}
	}
	else if(f->num_vertices() == 4){
		bool intersecting = false;
		if(RayTriangleIntersection(p, t1, t2, s, aaPos[f->vertex(0)], aaPos[f->vertex(1)],
									aaPos[f->vertex(2)], aaPos[e->vertex(0)], dir))
		{
			intersecting = true;
		}
		else if(RayTriangleIntersection(p, t1, t2, s, aaPos[f->vertex(0)], aaPos[f->vertex(2)],
										aaPos[f->vertex(3)], aaPos[e->vertex(0)], dir))
		{
			intersecting = true;
		}

		if(intersecting && (s >= 0) && (s <= 1.))
		{
		//	split the face
			vrt = *grid.create<RegularVertex>();
			Vertex* newEdgeVrts[4] = {NULL, NULL, NULL, NULL};
			refined = f->refine(newFaces, &newFaceVrt, newEdgeVrts, vrt);
		}
	}

	if(refined && vrt){
	//	create a new vertex and adjust position
		aaPos[vrt] = p;

	//	register the new faces and erase the old one
		for(size_t i = 0; i < newFaces.size(); ++i)
			grid.register_element(newFaces[i], f);
		grid.erase(f);

	//	to make sure that no double edges may occur, we'll use MergeVertices
	//	and SplitEdge
		RegularVertex* nVrt = SplitEdge<RegularVertex>(grid, grid, e);
		MergeVertices(grid, vrt, nVrt);

/*
	//	split the edge with the new vertex and erase the old one
		CreateEdgeSplitGeometry(grid, grid, e, vrt);
		grid.erase(e);
*/

		return true;
	}

	return false;
}

/**
 *	Projects vertices in elems onto close edges in elems.
 *	Though this method can be used to remove degenerated triangles,
 *	it is not guaranteed, that no degenerated triangles will remain
 *	(indeed, new degenerated triangles may be introduced).
 */
template <class TAAPosVRT>
bool ProjectVerticesToCloseEdges(Grid& grid,
								 GridObjectCollection elems,
								 TAAPosVRT& aaPos,
								 number snapThreshold)
{
//	perform vertex/edge intersections
//	iterate over all vertices
	for(VertexIterator vrtIter = elems.begin<Vertex>();
		vrtIter != elems.end<Vertex>();)
	{
		Vertex* vrt = *vrtIter;
		++vrtIter;

	//	check against all edges
		for(EdgeIterator eIter = elems.begin<Edge>();
			eIter != elems.end<Edge>();)
		{
			Edge* e = *eIter;
			++eIter;

		//	try to insert the vertex into the edge
			ResolveVertexEdgeIntersection(grid, vrt, e, aaPos, snapThreshold);
		}
	}

	return true;
}

/**
 *	Projects vertices in elems onto close faces in elems.
 */
template <class TAAPosVRT>
bool ProjectVerticesToCloseFaces(Grid& grid,
								 GridObjectCollection elems,
								 TAAPosVRT& aaPos,
								 number snapThreshold)
{
//	perform vertex/face intersections
//	iterate over all vertices
	for(VertexIterator vrtIter = elems.vertices_begin();
		vrtIter != elems.vertices_end();)
	{
		Vertex* vrt = *vrtIter;
		++vrtIter;

	//	check against all faces
		for(FaceIterator fIter = elems.faces_begin(); fIter != elems.faces_end();)
		{
			Face* f = *fIter;
			++fIter;

		//	try to insert the vertex into the face
			ResolveVertexFaceIntersection(grid, vrt, f, aaPos, snapThreshold);
		}
	}
	return true;
}

/**THIS METHOD USES Grid::mark.
 * Intersects all edges in elems which are closer to each other
 * than snapThreshold.*/
template <class TAAPosVRT>
bool IntersectCloseEdges(Grid& grid,
						 GridObjectCollection elems,
						 TAAPosVRT& aaPos,
						 number snapThreshold)
{
//	we'll first mark all elements in elems to make sure that
//	only edges which were initially marked are intersected.
	grid.begin_marking();
	for(EdgeIterator iter = elems.begin<Edge>();
		iter != elems.end<Edge>(); ++iter)
	{
		grid.mark(*iter);
	}

//	perform edge/edge and edge/face intersections
	for(EdgeIterator mainIter = elems.begin<Edge>();
		mainIter != elems.end<Edge>();)
	{
		Edge* e = *mainIter;
		++mainIter;

	//	if e is not marked, we can exit right away, since all succeeding
	//	edges won't be marked, too.
		if(!grid.is_marked(e))
			break;

	//	check all other edges up to e.
		for(EdgeIterator iter = elems.begin<Edge>(); *iter != e;)
		{
			Edge* e2 = *iter;
			++iter;

		//	if an intersection occured, we have to move on to the next edge in the queue,
		//	since the old edge no longer exists.
			if(ResolveEdgeEdgeIntersection(grid, e, e2, aaPos, snapThreshold)){
				break;
			}
		}
	}
	grid.end_marking();
	return true;
}


///	returns the index of the first vertex closer to p than snapThreshold.
/**	returns -1 if nothing was found.*/
template <class TAAPosVRT>
int FindCloseVertexInArray(std::vector<Vertex*>& array,
							const typename TAAPosVRT::ValueType& p,
							TAAPosVRT& aaPos, number snapThreshold)
{
	number snapThrSq = snapThreshold * snapThreshold;
//	iterate over the array and check whether a vertex close to vrt already exists.
	for(size_t i = 0; i < array.size(); ++i){
		if(VecDistanceSq(aaPos[array[i]], p) < snapThrSq){
		//	we got one. return the index
			return (int)i;
		}
	}
	return -1;
}

////////////////////////////////////////////////////////////////////////
/**	This method uses Grid::mark
 */
template <class TAAPosVRT>
bool ResolveGridIntersections(Grid& grid, TriangleIterator trisBegin,
							  TriangleIterator trisEnd, number snapThreshold,
							  TAAPosVRT& aaPos)
{
	using namespace std;
//todo: add octree
//	we use a selector to select elements that shall be merged and
//	triangles that are to be processed and deleted.
	Selector sel(grid);
	sel.enable_autoselection(false);

//	we first select all associated vertices and perform a merge on them
	sel.select(trisBegin, trisEnd);
	SelectAssociatedVertices(sel, trisBegin, trisEnd);
	SelectAssociatedEdges(sel, trisBegin, trisEnd);
	RemoveDoubles<3>(grid, sel.vertices_begin(), sel.vertices_end(),
					 aPosition, snapThreshold);

////////////////////////////////
//	PERFORM AND RESOLVE TRIANGLE - TRIANGLE INTERSECTIONS

//	clear edges and vertices from the selector. faces have to stay, since we will
//	operate on them now.
	sel.clear<Vertex>();
	sel.clear<Edge>();

//	enable selection inheritance, since we want new elements to be
//	selected in this selector
	sel.enable_selection_inheritance(true);

//	we need some attachments in order to store new vertices and edges for
//	each face.
	typedef Attachment<vector<Vertex*> >		AVrtVec;
	typedef Attachment<vector<pair<int, int> > >	AEdgeDescVec;
	AVrtVec aVrtVec;
	AEdgeDescVec aEdgeDescVec;
	grid.attach_to_faces(aVrtVec);
	grid.attach_to_faces(aEdgeDescVec);
	Grid::FaceAttachmentAccessor<AVrtVec> aaVrtVec(grid, aVrtVec);
	Grid::FaceAttachmentAccessor<AEdgeDescVec> aaEdgeDescVec(grid, aEdgeDescVec);

//	iterate over all triangles and perform intersecion with other triangles
	for(TriangleIterator triIter1 = sel.begin<Triangle>();
		triIter1 != sel.end<Triangle>(); ++triIter1)
	{
		Triangle* t1 = *triIter1;

	//	iterate over the rest of the triangles
		TriangleIterator triIter2 = triIter1;
		for(++triIter2; triIter2 != sel.end<Triangle>(); ++triIter2)
		{
			Triangle* t2 = *triIter2;

		//	we have to make sure, that t1 and t2 do not share an edge (two vertices)
			size_t numShared = NumSharedVertices(grid, t1, t2);
			if(numShared > 1)
				continue;

		//	perform normal comparision to avoid intersection of flat neighbours
			vector3 n1, n2;
			CalculateNormal(n1, t1, aaPos);
			CalculateNormal(n2, t2, aaPos);
			number d = VecDot(n1, n2);
			if(fabs(d) > 1. - SMALL)
				continue;

			vector3 ip[2];
			if(TriangleTriangleIntersection(aaPos[t1->vertex(0)], aaPos[t1->vertex(1)],
											aaPos[t1->vertex(2)], aaPos[t2->vertex(0)],
											aaPos[t2->vertex(1)], aaPos[t2->vertex(2)],
											&ip[0], &ip[1]) == 1)
			{
			//	add an edge between the two points
			//	to avoid insertion of double points, we first check whether the point
			//	already exists in the triangle. Do this for both triangles.
				Triangle* t[2]; t[0] = t1; t[1] = t2;

			//	prepare both triangles.
			//todo: think about performance optimizations.
			//	insertion of corner points could be avoided by bloating the code a little.
			//	this could increase performance.
				for(size_t i_tri = 0; i_tri < 2; ++i_tri){
				//	If it is encountered for the first time,
				//	we'll add its corner-vertices to its list of vertices.
					vector<Vertex*>& vrts = aaVrtVec[t[i_tri]];
					if(vrts.empty()){
						for(size_t i = 0; i < t[i_tri]->num_vertices(); ++i)
							vrts.push_back(t[i_tri]->vertex(i));
					}
				}

				//	now check whether the vertex already exists
				int inds1[2];
				int inds2[2];
				for(size_t i = 0; i < 2; ++i){
					int tind1 = FindCloseVertexInArray(aaVrtVec[t[0]], ip[i],
													   aaPos, snapThreshold);
					int tind2 = FindCloseVertexInArray(aaVrtVec[t[1]], ip[i],
													   aaPos, snapThreshold);

					if(tind1 == -1){
						if(tind2 == -1){
						//	we have to create a new vertex
							Vertex* vrt = *grid.create<RegularVertex>();
							aaPos[vrt] = ip[i];
							tind1 = (int)aaVrtVec[t[0]].size();
							tind2 = (int)aaVrtVec[t[1]].size();
							aaVrtVec[t[0]].push_back(vrt);
							aaVrtVec[t[1]].push_back(vrt);
						}
						else{
						//	the vertex already exists in t[1]
							tind1 = (int)aaVrtVec[t[0]].size();
							aaVrtVec[t[0]].push_back((aaVrtVec[t[1]])[tind2]);
						}
					}
					else if(tind2 == -1){
					//	the vertex already exists in t[0]
						tind2 = (int)aaVrtVec[t[1]].size();
						aaVrtVec[t[1]].push_back((aaVrtVec[t[0]])[tind1]);
					}

				//	ind1 now contains the index into the vertex array of t[0], at
				//	which a vertex with position ip[i] lies.
					inds1[i] = tind1;
					inds2[i] = tind2;
				}

			//	we found the indices of both endpoints and can now add an edge
			//	connecting both to the edgeDesc arrays of t[0] and t[1].
				if(inds1[0] != inds1[1])
					aaEdgeDescVec[t[0]].push_back(make_pair(inds1[0], inds1[1]));
				if(inds2[0] != inds2[1])
					aaEdgeDescVec[t[1]].push_back(make_pair(inds2[0], inds2[1]));
			}
		}
	}

//	all intersections have been resolved. Iterate over the triangles again and
//	create the new elements.
//	triangles that shall be deleted are pushed to vDelTris
	vector<Triangle*> vDelTris;
//	here we collect all vertices on which a merge has to be performed at the end
//	of the algorithm (vertices created through edge-edge intersections inside a triangle)
	vector<Vertex*> cutVertices;
	Grid tgrid(GRIDOPT_STANDARD_INTERCONNECTION);
	AInt aInt;
	AVertex aVrt;
	tgrid.attach_to_vertices(aPosition);
	tgrid.attach_to_vertices(aInt);
	tgrid.attach_to_vertices_dv(aVrt, NULL);
	Grid::VertexAttachmentAccessor<APosition> taaPos(tgrid, aPosition);
	Grid::VertexAttachmentAccessor<AVertex> aaVrt(tgrid, aVrt);

//	holds vertices of tgrid, so that they are accessible by index.
	vector<Vertex*> tgridVrts;

	for(TriangleIterator triIter = sel.begin<Triangle>();
		triIter != sel.end<Triangle>(); ++triIter)
	{
		Triangle* t = *triIter;

	//	we only proceed if there are intersecion-edges at all
		if(!aaEdgeDescVec[t].empty()){
			tgrid.clear_geometry();
			tgridVrts.clear();

		//	copy vertices associated with t1 to tgrid
			vector<Vertex*>& vrts = aaVrtVec[t];
			for(size_t i = 0; i < vrts.size(); ++i){
				Vertex* vrt = *tgrid.create<RegularVertex>();
				aaVrt[vrt] = vrts[i];
				taaPos[vrt] = aaPos[vrts[i]];
				tgridVrts.push_back(vrt);
			}

		//	now create the edges. vertices are found by indexing tgridVrts
			vector<pair<int, int> >& edgeDescs = aaEdgeDescVec[t];

		//	tri edges
			tgrid.create<RegularEdge>(EdgeDescriptor(tgridVrts[0], tgridVrts[1]));
			tgrid.create<RegularEdge>(EdgeDescriptor(tgridVrts[1], tgridVrts[2]));
			tgrid.create<RegularEdge>(EdgeDescriptor(tgridVrts[2], tgridVrts[0]));

		//	new edges
			for(size_t i = 0; i < edgeDescs.size(); ++i){
				tgrid.create<RegularEdge>(EdgeDescriptor(tgridVrts[edgeDescs[i].first],
												  tgridVrts[edgeDescs[i].second]));
			}

		//	we now have to resolve intersections between the edges
		//	first we'll try to snap vertices to edges
			ProjectVerticesToCloseEdges(tgrid, tgrid.get_grid_objects(),
										taaPos, SMALL);

		//	now resolve edge/edge intersections
			IntersectCloseEdges(tgrid, tgrid.get_grid_objects(),
								taaPos, SMALL);

		//	make sure that all vertices have an associated aaVrt
			for(VertexIterator viter = tgrid.vertices_begin();
				viter != tgrid.vertices_end(); ++viter)
			{
				if(!aaVrt[*viter]){
				//	since the vertex does not have an associated vertex in grid,
				//	it is clear that it has been created through an edge-edge cut.
				//	Associates of such vertices have to be merged later on.
					aaVrt[*viter] = *grid.create<RegularVertex>();
					aaPos[aaVrt[*viter]] = taaPos[*viter];
					cutVertices.push_back(aaVrt[*viter]);
				}
			}

		//	ok. Everything is prepared. We can now triangulate the grid.
			if(TriangleFill_SweepLine(tgrid, tgrid.edges_begin(), tgrid.edges_end(),
										aPosition, aInt))
			{
			//	mark the triangle for deletion
				vDelTris.push_back(*triIter);

			//	add the triangles to the grid.
				for(TriangleIterator titer = tgrid.begin<Triangle>();
					titer != tgrid.end<Triangle>(); ++titer)
				{
					Triangle* ntri = *titer;

					grid.create<Triangle>(TriangleDescriptor(aaVrt[ntri->vertex(0)],
															aaVrt[ntri->vertex(1)],
															aaVrt[ntri->vertex(2)]),
										 *triIter);
				}
			}
			else{/*
				static int fileCounter = 1;
				string filenamePrefix = "/Users/sreiter/Desktop/failed_sweeplines/failed_sweepline_";
				stringstream ss2d, ss3d;
				ss2d << filenamePrefix << "2d_" << fileCounter << ".lgb";
				ss3d << filenamePrefix << "3d_" << fileCounter << ".lgb";
				++fileCounter;
				UG_LOG("TriangleFill_SweepLine failed!\n");
				SaveGridToFile(tgrid, ss3d.str().c_str(), aPosition);
			//	perform transformation to 2d and save that too.
				std::vector<vector3> vrts;
				for(VertexIterator iter = tgrid.vertices_begin();
					iter != tgrid.vertices_end(); ++iter)
				{
					vrts.push_back(taaPos[*iter]);
				}
				std::vector<vector2> vrts2d(vrts.size());
				TransformPointSetTo2D(&vrts2d.front(), &vrts.front(),
									  vrts.size());

				size_t counter = 0;
				for(VertexIterator iter = tgrid.vertices_begin();
					iter != tgrid.vertices_end(); ++iter, ++counter)
				{
					taaPos[*iter] = vector3(vrts2d[counter].x(), vrts2d[counter].y(), 0);
				}

				SaveGridToFile(tgrid, ss2d.str().c_str(), aPosition);
				*/
			}
		}
	}

//	detach attachments (tgrid is deleted anyways)
	grid.detach_from_faces(aVrtVec);
	grid.detach_from_faces(aEdgeDescVec);

////////////////////////////////
//	GRID POSTPROCESS
//	before we merge vertices in cutVertices, we'll select all faces
//	in order to make sure that only valid faces will be deleted.
	sel.clear();
	sel.select(vDelTris.begin(), vDelTris.end());
	sel.select(cutVertices.begin(), cutVertices.end());

//	perform the merge (this has to be done on a selector.
//	  the current version of RemoveDoubles is a little restrictive
//	  in this regard.)
	if(!sel.empty<Vertex>()){
		RemoveDoubles<3>(grid, sel.vertices_begin(), sel.vertices_end(),
						 aPosition, snapThreshold);
	}

	sel.clear<Vertex>();
	sel.clear<Edge>();

//	finally delete all refined triangles and associated unused edges and vertices
	SelectInnerSelectionEdges(sel);
	SelectInnerSelectionVertices(sel);

	grid.erase(sel.begin<Face>(), sel.end<Face>());
	grid.erase(sel.begin<Edge>(), sel.end<Edge>());
	grid.erase(sel.begin<Vertex>(), sel.end<Vertex>());

	return true;
}

}// end of namespace

#endif
