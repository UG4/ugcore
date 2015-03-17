// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include "delaunay_triangulation.h"
#include "lib_grid/algorithms/grid_generation/triangle_fill_sweep_line.h"
#include "lib_grid/algorithms/subset_util.h"
#include "lib_grid/algorithms/polychain_util.h"
//temporary
#include "lib_grid/file_io/file_io.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
class DelaunayDebugSaver
{
	public:
		static DelaunayDebugSaver& inst()
		{
			static DelaunayDebugSaver dds;
			return dds;
		}

		void save(Grid& g, const char* msg)
		{
			if(m_saveEnabled){
				std::stringstream ss;
				ss << "dbg" << m_counter++ << ".ugx";

				UG_LOG(msg << ": debug-save to " << ss.str() << std::endl);
				SubsetHandler sh(g);
				AssignGridToSubset(g, sh, 0);
				SaveGridToFile(g, sh, ss.str().c_str());
			}
		}

		template <class TAAPos>
		void save(Grid& g, const char* msg, DelaunayInfo<TAAPos>& dinfo)
		{
			if(m_saveEnabled){
				std::stringstream ss;
				ss << "dbg" << m_counter++ << ".ugx";

				UG_LOG(msg << ": debug-save to " << ss.str() << std::endl);
				SubsetHandler sh(g);

				for(VertexIterator iter = g.begin<Vertex>();
					iter != g.end<Vertex>(); ++iter)
				{
					if(dinfo.is_marked(*iter))
						sh.assign_subset(*iter, 1);
					else
						sh.assign_subset(*iter, 0);
				}

				for(EdgeIterator iter = g.begin<Edge>();
					iter != g.end<Edge>(); ++iter)
				{
					if(dinfo.is_constrained(*iter))
						sh.assign_subset(*iter, 3);
					if(dinfo.is_candidate(*iter))
						sh.assign_subset(*iter, 4);
					else
						sh.assign_subset(*iter, 2);
				}

				for(FaceIterator iter = g.begin<Face>();
					iter != g.end<Face>(); ++iter)
				{
					if(dinfo.is_marked(*iter))
						sh.assign_subset(*iter, 6);
					else
						sh.assign_subset(*iter, 5);
				}

				sh.subset_info(0).name = "unmarked vrts";
				sh.subset_info(1).name = "marked vrts";
				sh.subset_info(2).name = "unmarked edges";
				sh.subset_info(3).name = "constrained edges";
				sh.subset_info(4).name = "candidate edges";
				sh.subset_info(5).name = "unmarked faces";
				sh.subset_info(6).name = "marked faces";

				AssignSubsetColors(sh);
				SaveGridToFile(g, sh, ss.str().c_str());

			}
		}

		void enable_save(bool enable = true)	{m_saveEnabled = enable;}

	private:
		DelaunayDebugSaver() : m_saveEnabled(false), m_counter(0) {}

		bool	m_saveEnabled;
		int		m_counter;
};

static void EnableDelaunayDebugSave(bool enable = true)
{
	DelaunayDebugSaver::inst().enable_save(enable);
}

// static void DelaunayDebugSave(Grid& g, const char* msg)
// {
// 	DelaunayDebugSaver::inst().save(g, msg);
// }

template <class TAAPos>
static void DelaunayDebugSave(Grid& g, const char* msg, DelaunayInfo<TAAPos>& dinfo)
{
	DelaunayDebugSaver::inst().save(g, msg, dinfo);
}



////////////////////////////////////////////////////////////////////////////////
template <class TAAPos>
bool MakeDelaunay(DelaunayInfo<TAAPos>& info)
{
	using namespace std;
	typedef typename TAAPos::ValueType vector_t;

	Grid& grid = info.grid();
	Face* nbrFaces[2];
	vector<Edge*> edges;

//UG_LOG("MakeDelaunay...\n");
	TAAPos& aaPos = info.position_accessor();

	while(info.candidates_left()){
		Edge* e = info.pop_candidate();

	//	we only perform swaps on regular manifolds.
		if(GetAssociatedFaces(nbrFaces, grid, e, 2) == 2){
		//	make sure that both neighbors are triangles
			if(nbrFaces[0]->num_vertices() != 3 || nbrFaces[1]->num_vertices() != 3)
				continue;

		//	This section is just temporary...
			Vertex* conVrt0 = GetConnectedVertex(e, nbrFaces[0]);
			Vertex* conVrt1 = GetConnectedVertex(e, nbrFaces[1]);

			vector_t& v0 = aaPos[e->vertex(0)];
			vector_t& v1 = aaPos[e->vertex(1)];
			vector_t& v2 = aaPos[conVrt0];
			vector_t& v3 = aaPos[conVrt1];

			vector_t cc1, cc2;

			bool cc1_ok = TriangleCircumcenter(cc1, v0, v1, v2);
			bool cc2_ok = TriangleCircumcenter(cc2, v1, v0, v3);

			number r1Sq = VecDistanceSq(cc1, v0);
			number r2Sq = VecDistanceSq(cc2, v0);

//			UG_LOG("  checking edge: " << CalculateCenter(e, aaPos) << "\n");
/*					<< "cc1_ok: " << cc1_ok << ", cc2_ok: " << cc2_ok
					<< ", cc1: " << cc1 << ", cc2: " << cc2
					<< ", r1: " << sqrt(r1Sq) << ", r2: " << sqrt(r2Sq) << endl);
*/

		//	for stability reasons, we're checking against the smaller circle
			if(cc1_ok){
				if(cc2_ok){
					if(r1Sq <= r2Sq){
						if(r1Sq <= VecDistanceSq(cc1, v3))
							continue; //the edge is fine
					}
					else{
						if(r2Sq <= VecDistanceSq(cc2, v2))
							continue; //the edge is fine
					}
				}
				else{
					if(r1Sq <= VecDistanceSq(cc1, v3))
						continue; //the edge is fine
				}
			}
			else if(cc2_ok){
				if(r2Sq <= VecDistanceSq(cc2, v2))
					continue; //the edge is fine
			}
			else{
				UG_LOG("TriangleCircumcenter failed! Excpect non-delaunay output!\n");
				UG_LOG("  This is most likely caused by two degenerated triangles which "
						"share an edge.\n");
				UG_LOG("edge center: " << CalculateCenter(e, aaPos) << endl);
				return false;
				//continue;
			}

//			UG_LOG("  swap-it!\n");

		//	before swapping, we have to make sure, that the generated triangle
		//	won't be degenerated.
		//	this is a costly operation... we check whether both circumcenters
		//	of the new triangles can be calculated...
/*
			if(!(TriangleCircumcenter(cc1, v0, v2, v3)
				 && TriangleCircumcenter(cc2, v2, v1, v3)))
			{
			//	we have to abort the swap
				continue;
			}
*/
		//	ok - everything is fine. Now swap the edge
			Edge* eNew = SwapEdge(grid,  e);

			if(!eNew){
				UG_LOG("An edge-swap failed. Expect degenerated or flipped triangles "
						"and a non-delaunay output!\n");
				UG_LOG("edge center: " << CalculateCenter(e, aaPos) << endl);
				return false;
				//continue;
			}

			e = eNew;

			//DelaunayDebugSave(grid, "Edge Swapped", info);

		//	all edges of associated triangles are candidates again (except e)
			GetAssociatedFaces(nbrFaces, grid, e, 2);
			for(size_t i = 0; i < 2; ++i){
				CollectAssociated(edges, grid, nbrFaces[i]);
				for(size_t j = 0; j < edges.size(); ++j){
					if(edges[j] != e && !(info.is_constrained(edges[j])
										  || info.is_candidate(edges[j])))
					{
						info.push_candidate(edges[j]);
					}
				}
			}
		}
	}
	return true;
}


static bool
DelaunayLineLineIntersection(vector3& vOut,
							 const vector3& lineFrom, const vector3& lineTo,
							 const vector3& edgeVrt1, const vector3& edgeVrt2,
							 vector3 areaNormal)
{
	number t;

	vector3 lineDir, edgeDir, planeNormal;
	VecSubtract(lineDir, lineTo, lineFrom);
	VecNormalize(lineDir, lineDir);
	VecSubtract(edgeDir, edgeVrt2, edgeVrt1);
	VecNormalize(edgeDir, edgeDir);
	VecCross(planeNormal, edgeDir, areaNormal);
	if(RayPlaneIntersection(vOut, t, lineFrom, lineDir, edgeVrt1, planeNormal)){
		if(t > 0 && VecDistanceSq(lineFrom, vOut) < VecDistanceSq(lineFrom, lineTo) + SMALL_SQ){
			return true;
		}
	}
	return false;
}



template <class TAAPos>
bool QualityGridGeneration(Grid& grid, DelaunayInfo<TAAPos>& info,
						   number minAngle, int maxSteps)
{
	EnableDelaunayDebugSave(false);
	// maxSteps = 100;
	// UG_LOG("DEBUG: Setting maxSteps to " << maxSteps << "\n");

	using namespace std;
	typedef typename TAAPos::ValueType vector_t;

	TAAPos aaPos = info.position_accessor();

	int stepCount = 0;

//	helper to collect neighbors
	Face* nbrFaces[2];
	queue<Vertex*> qvrts; // used during splits
	vector<Vertex*> vrts;
	vector<Edge*> edges;
	vector<Edge*> closeEdges; // used during splits
	vector<Face*> faces;

	//UG_LOG("MakeDelaunay Initial\n");
	MakeDelaunay(info);

	if(minAngle > 0 && maxSteps != 0){
		info.enable_face_classification(minAngle);

	//	while there are faces left which have to be improved
		while(info.classified_faces_left()){
			/*if(stepCount == 488){
				EnableDelaunayDebugSave();
			}*/

			Face* f = info.pop_classified_face();
			if(f->num_vertices() != 3)
				continue;

		//todo: the normal is only required for 3d-types. indeed this will crash
		//		for 2d. This should be moved to a Ray-Line_Intersection3d test.
			vector3 faceNormal;
			CalculateNormal(faceNormal, f, aaPos);
		//	we can't operate on degenerated faces. Let's hope, that this face
		//	will improve during improvement of some non-degenerated face.
			if(VecLengthSq(faceNormal) < SMALL)
				return false;

			vector_t faceCenter = CalculateCenter(f, aaPos);

		//	if the triangle has 3 constrained edges, we'll ignore it.
		//	if it has two, we'll ignore it if the angle between the
		//	two constrained ones is smaller than minAngle.
			{
				CollectAssociated(edges, grid, f);
				int numConstrained = 0;
				Edge* ce[3];
				for(size_t i = 0; i < edges.size(); ++i){
					if(info.is_constrained(edges[i])){
						ce[numConstrained] = edges[i];
						++numConstrained;
					}
				}
				if(numConstrained > 1)
					continue;
				else if(numConstrained == 2){
					// UG_LOG("Examining face with two constrained edges at: " << faceCenter << endl);
					Vertex* sharedVrt = NULL, *v1 = NULL, *v2 = NULL;
					for(int i = 0; (i < 2) && (!sharedVrt); ++i){
						for(int j = 0; j < 2; ++j){
							if(ce[0]->vertex(i) == ce[1]->vertex(j)){
								sharedVrt = ce[0]->vertex(i);
								v1 = ce[0]->vertex((i+1)%2);
								v2 = ce[1]->vertex((j+1)%2);
								break;
							}
						}
					}
					if(!sharedVrt)
						continue;
					vector_t d1, d2;
					VecSubtract(d1, aaPos[v1], aaPos[sharedVrt]);
					VecSubtract(d2, aaPos[v2], aaPos[sharedVrt]);
					if(rad_to_deg(VecAngle(d1, d2)) <= minAngle)
						continue;
				}
			}

		//	calculate triangle-circumcenter
			vector_t& v0 = aaPos[f->vertex(0)];
			vector_t& v1 = aaPos[f->vertex(1)];
			vector_t& v2 = aaPos[f->vertex(2)];
			vector_t cc;

			if(!TriangleCircumcenter(cc, v0, v1, v2)){
				UG_LOG("Couldn't calculate triangle-circumcenter. Expect unexpected results!\n");
				UG_LOG("triangle: " << faceCenter << "\n");
				//SaveGridToFile(grid, "delaunay_debug.ugx");
				return false;
				//continue;
			}

		//	locate the triangle which contains cc. Do this by traversing edges
		//	as required. Note that since the delaunay property holds, we're only
		//	traversing edges in a circle, which does not contain any vertices.
			Edge* lastTraversedEdge = NULL;
			Face* curFace = f;
			//vector_t startPos = CalculateCenter(f, aaPos);
			vector_t rayDir;
			VecSubtract(rayDir, cc, faceCenter);
			Vertex* pointInserted = NULL;

			// UG_LOG("curFace: " << faceCenter << "\n");
			while(pointInserted == NULL){
				//UG_LOG("curTri: " << CalculateCenter(curFace, aaPos) << "\n");
			//	to make things as robust as possible, we'll always intersect the
			//	same line with each edge
				Edge* nextEdge = NULL;
				bool split = false;

				CollectAssociated(edges, grid, curFace);
				for(size_t i = 0; i < edges.size(); ++i){
					Edge* e = edges[i];
					if(e == lastTraversedEdge)
						continue;

					vector_t a;
					if(DelaunayLineLineIntersection(a, faceCenter, cc,
							aaPos[e->vertex(0)], aaPos[e->vertex(1)], faceNormal))
					{
					//	this edge will be traversed next.
						nextEdge = e;
					//	check whether e has to be split, to avoid bad aspect ratios
						split = (VecDistanceSq(faceCenter, a) > VecDistanceSq(faceCenter, cc) - SMALL_SQ);
						break;
					}/*
					// This else-condition is not necessary - if they are parallel, that's not a problem...
					else{
						UG_LOG("Ray-Plane intersection failed in step " << stepCount << "... aborting.\n");
						UG_LOG("  curTri: " << CalculateCenter(curFace, aaPos)
								<< ", curEdge: " << CalculateCenter(e, aaPos) << endl);
						UG_LOG("  face-normal: " << faceNormal << ", plane-normal: " << planeNormal
								<< "  ray-dir: " << rayDir << endl);
						return false;
					}*/
				}

				if(nextEdge){
				//	check whether the edge has to be splitted
					bool isConstrained = info.is_constrained(nextEdge);
					split |= isConstrained;
					if(!split){
						int numNbrs = GetAssociatedFaces(nbrFaces, grid, nextEdge, 2);
						if(numNbrs != 2)
							split = true;
						else{
						//	get the next face
							if(nbrFaces[0] == curFace)
								curFace = nbrFaces[1];
							else
								curFace = nbrFaces[0];
						//	and set the last traversed edge
							lastTraversedEdge = nextEdge;
						}
					}

					if(split){
						//UG_LOG("EDGE-SPLIT\n");
						vector_t center = CalculateCenter(nextEdge, aaPos);

						Vertex* vrt0 = nextEdge->vertex(0);
						Vertex* vrt1 = nextEdge->vertex(1);
						number radiusSq = VecDistanceSq(center, aaPos[vrt0]);
						number radius = sqrt(radiusSq);
						pointInserted = SplitEdge<RegularVertex>(grid, nextEdge, false);

						if(isConstrained){
							// UG_LOG("IS CONSTRAINED\n");
							aaPos[pointInserted] = center;

						//	we have to erase all vertices, which are marked and
						//	which are closer to the inserted point than the
						//	radius calculated above.
						//todo: And which are visible (not hidden by constrained edges)
							assert(qvrts.empty());
							qvrts.push(pointInserted);

							vrts.clear();
							closeEdges.clear();

							grid.begin_marking();
							while(!qvrts.empty()){
								Vertex* vrt = qvrts.front();
								qvrts.pop();

							//	collect all edges connected to vrt
								CollectAssociated(edges, grid, vrt);
							//	if the edge is not yet marked, we'll add it to the
							//	closeEdges list.
								for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
									Edge* e = edges[i_edge];
									if(!grid.is_marked(e)){
									//	check whether the edge intersects the critical circle
										if(DistancePointToLine(center, aaPos[e->vertex(0)], aaPos[e->vertex(1)])
											< radius)
										{
										//	if the edge is a constrained edge, we'll push it to
										//	closeEdges
											if(info.is_constrained(e) && !EdgeContains(e, pointInserted))
												closeEdges.push_back(e);

											grid.mark(e);

										//	the connected vertex has to be pushed to
										//	the queue, regardless whether it lies in the circle
										//	or not, since an edge from this vertex could
										//	reenter the circle
											Vertex* vcon = GetConnectedVertex(e, vrt);
											if(grid.is_marked(vcon))
												continue;
											grid.mark(vcon);
											qvrts.push(vcon);

										//	vrt0 and vrt1 won't be removed anyways
											if(vcon == vrt0 || vcon == vrt1)
												continue;

										//	if the new vertrex is marked (created during remeshing)
										//	and lies in the circle, then we'll push it to vrts
											//UG_LOG("is_marked: " << info.is_marked(vcon) << "rad: " << radiusSq
											//	<< ", testvrt: " << VecDistanceSq(aaPos[vcon], center) << "\n");
											if(VecDistanceSq(aaPos[vcon], center) < radiusSq){
												if(info.is_marked(vcon))
													vrts.push_back(vcon);
											}
										}
									}
								}
							}
							grid.end_marking();

							//UG_LOG(vrts.size() << " vertices have to be deleted!\n");

						//	vrts now contains all vertices which lie in a circle of radius
						//	'radiusSq' around the center. We now have to check for each,
						//	whether it is visible from the center, by checking whether
						//	the line of sight intersects a constrained edge in
						//	closeEdges.
						//	if it is visible, then we'll delete it and locally remesh
						//	the grid.
						//	We want to record all new edges that are created during remeshing,
						//	since they have to be considered as candidates for MakeDelaunay.
						//	However, since some of them will possibly be deleted during
						//	following erasures, we have to be careful here. This is
						//	all handled by info.start/stop_candidate_recording
							info.start_candidate_recording();
							for(size_t i_vrts = 0; i_vrts < vrts.size(); ++i_vrts){
								Vertex* vrt = vrts[i_vrts];
								vector_t& v = aaPos[vrt];
								//UG_LOG("ERASING VRT at " << v << "...\n");
								bool intersects = false;
								//UG_LOG("  checking for intersections\n");
								for(size_t i_edge = 0; i_edge < closeEdges.size(); ++i_edge){
									vector_t& ev0 = aaPos[closeEdges[i_edge]->vertex(0)];
									vector_t& ev1 = aaPos[closeEdges[i_edge]->vertex(1)];
									vector_t a;
									if(DelaunayLineLineIntersection(a, faceCenter, v, ev0, ev1, faceNormal)){
									//if(LineLineIntersection3d(a, b, center, v, ev0, ev1)){
										//UG_LOG("  intersection!\n");
										intersects = true;
										break;
									}
								}

								if(intersects)
									continue;

								//UG_LOG("  gathering surrounding edges\n");
							//	no intersection detected. erase the vertex and perform
							//	local retriangulation. Add surrounding edges to
							//	the 'edges' array
							//	Store one face, which will be the parent for all new faces
								Face* parent = NULL;
								edges.clear();
								CollectAssociated(faces, grid, vrt);
								for(size_t i_face = 0; i_face < faces.size(); ++i_face){
									Face* f = faces[i_face];
									parent = f;
									for(size_t i = 0; i < f->num_edges(); ++i){
										Edge* e = grid.get_edge(f, i);
										if(!EdgeContains(e, vrt)){
											edges.push_back(e);
											//UG_LOG("surrounding edge: " << CalculateCenter(e, aaPos) << endl);
										}
									}
								}
								assert(parent);

								//UG_LOG("  retriangulation\n");
							//	perform retriangulation
							//	get the vertices of the poly-chain
								//UG_LOG("1");
								std::vector<Vertex*> vrtPolyChain;
								CreatePolyChain(vrtPolyChain, grid, edges.begin(), edges.end());
								//UG_LOG("2");
								std::vector<vector_t> posPolyChain(vrtPolyChain.size());
								std::vector<int> edgeChain, faceIndices;
								edgeChain.reserve(vrtPolyChain.size() * 2);
								//UG_LOG("3");
								size_t numVrts = vrtPolyChain.size();
								for(size_t i = 0; i < numVrts; ++i){
									//UG_LOG(aaPos[vrtPolyChain[i]]);
									posPolyChain[i] = aaPos[vrtPolyChain[i]];
									edgeChain.push_back(i);
									edgeChain.push_back((i+1) % numVrts);
								}
/*
								("edges.size: " << edges.size() << ", posPolyChain.size: " << posPolyChain.size());
								for(size_t i = 0; i < edgeChain.size(); i+=2){
									UG_LOG(" (" << edgeChain[i] << ", " << edgeChain[i+1] << ")");
								}
								UG_LOG(endl);
*/

								//UG_LOG("4");

								if(!TriangleFill_SweepLine(faceIndices, posPolyChain, edgeChain)){
									UG_LOG("Sweepline failed in step " << stepCount << "... aborting\n");
									UG_LOG("  While examining face " << faceCenter << endl);
									return false;
								}
								//UG_LOG("5");
							//	create the faces
								for(size_t i = 0; i < faceIndices.size(); i+=3){
									grid.create<Triangle>(
											TriangleDescriptor(vrtPolyChain[faceIndices[i]],
															   vrtPolyChain[faceIndices[i+1]],
															   vrtPolyChain[faceIndices[i+2]]),
											parent);
								}
								//UG_LOG("6\n");

								//UG_LOG("  erase\n");
							//	now erase vrt
								grid.erase(vrt);

								//DelaunayDebugSave(grid, "After TriangleFill", info);
							}

							info.stop_candidate_recording();
						}
						else{
							aaPos[pointInserted] = cc;
						}
					}
				}
				else{
				//	we found the triangle, which contains cc. Insert the point
				//	and perform local delaunay (not necessarily local...).
				//	todo: make sure, that cc really lies in curTri
				//...
					//UG_LOG("Inserting point into triangle\n");
					Vertex* vrt0 = curFace->vertex(0);
					Vertex* vrt1 = curFace->vertex(1);
					Vertex* vrt2 = curFace->vertex(2);

				//UG_LOG("creating new elements...\n");
					Vertex* vrt = *grid.create<RegularVertex>(curFace);
					aaPos[vrt] = cc;
					grid.create<Triangle>(TriangleDescriptor(vrt0, vrt1, vrt), curFace);
					grid.create<Triangle>(TriangleDescriptor(vrt1, vrt2, vrt), curFace);
					grid.create<Triangle>(TriangleDescriptor(vrt2, vrt0, vrt), curFace);
				//UG_LOG("erasing old face\n");
					grid.erase(curFace);

					//DelaunayDebugSave(grid, "After InsertPoint", info);

					pointInserted = vrt;
				}

				if(pointInserted){
					//UG_LOG("temp-save to delaunay_debug.ugx\n");
					//SaveGridToFile(grid, "delaunay_debug.ugx");

					//UG_LOG("inserted point\n");
				//	if a vertex was inserted, we'll have to perform a delaunay step.
				//	Find new candidates by examining edges of associated triangles of vrt.
				//	If an edge of such a triangle is connected to exactly 2
				//	triangles, both marked, then it is a new candidate.
					CollectAssociated(faces, grid, pointInserted);
					for(size_t i_face = 0; i_face < faces.size(); ++i_face){
						if(!info.is_marked(faces[i_face]))
							continue;

						CollectAssociated(edges, grid, faces[i_face]);
						for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
							Edge* e = edges[i_edge];
							if(info.is_candidate(e) || info.is_constrained(e))
								continue;

							Face* nbrs[2];
							if(GetAssociatedFaces(nbrs, grid, e, 2) == 2){
								if(info.is_marked(nbrs[0])
								   && info.is_marked(nbrs[1]))
								{
									info.push_candidate(e);
								}
							}
						}
					}

					// DelaunayDebugSave(grid, "Candidates Adjusted After Insert Point", info);

					//UG_LOG("redelaunaylizing\n");
					// UG_LOG("MakeDelaunay after InsertPoint\n");
					if(!MakeDelaunay(info)){
						UG_LOG("Make Delaunay failed in step " << stepCount << ".\n");
						UG_LOG("  While examining face " << faceCenter << endl);
						return false;
					}
					// UG_LOG("MakeDelaunay done\n");
					// UG_LOG("adaption step completed!\n");
				}
			}

/*
		//	check whether an illegal triangle has been inserted
			for(TriangleIterator iter = grid.begin<Triangle>();
				iter != grid.end<Triangle>(); ++iter)
			{
				vector3 n;
				CalculateNormal(n, *iter, aaPos);
				if(n.z() <= 0){
					UG_LOG("ATTENTION: Illegal triangle created in step " << stepCount << "\n");
					return false;
				}
			}
*/
			++stepCount;
			DelaunayDebugSave(grid, "", info);
			if(stepCount == maxSteps)
				break;
/*
			if((stepCount % 1000) == 0){
				UG_LOG("temp-save to delaunay_debug.ugx in step " << stepCount << endl);
				//DelaunayDebugSave(grid, "", info);
				SaveGridToFile(grid, "delaunay_debug.ugx");
			}
*/
		}
	}
	return true;
}



template bool MakeDelaunay(DelaunayInfo<Grid::VertexAttachmentAccessor<AVector3> >&);

template bool QualityGridGeneration(Grid&,
						   DelaunayInfo<Grid::VertexAttachmentAccessor<AVector3> >&,
						   number, int);

}//	end of namespace
