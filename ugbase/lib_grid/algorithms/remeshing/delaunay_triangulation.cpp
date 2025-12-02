/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#include <algorithm>
#include "delaunay_triangulation.h"
#include "lib_grid/algorithms/grid_generation/triangle_fill_sweep_line.h"
#include "lib_grid/algorithms/subset_util.h"
#include "lib_grid/algorithms/polychain_util.h"
//temporary
#include "lib_grid/file_io/file_io.h"

using namespace std;

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

		template <typename TAAPos>
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
					sh.assign_subset(*iter, dinfo.mark(*iter));
				}

				for(EdgeIterator iter = g.begin<Edge>();
					iter != g.end<Edge>(); ++iter)
				{
					sh.assign_subset(*iter, dinfo.mark(*iter));
				}

				for(FaceIterator iter = g.begin<Face>();
					iter != g.end<Face>(); ++iter)
				{
					sh.assign_subset(*iter, dinfo.mark(*iter));
				}

				sh.subset_info(0).name = "none";
				sh.subset_info(1).name = "inner";
				sh.subset_info(2).name = "new-inner";
				sh.subset_info(3).name = "segment";
				sh.subset_info(4).name = "new-segment";
				sh.subset_info(5).name = "dart";
				sh.subset_info(6).name = "shell";

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

template <typename TAAPos>
static void DelaunayDebugSave(Grid& g, const char* msg, DelaunayInfo<TAAPos>& dinfo)
{
	DelaunayDebugSaver::inst().save(g, msg, dinfo);
}



////////////////////////////////////////////////////////////////////////////////
template <typename TAAPos>
bool MakeDelaunay(DelaunayInfo<TAAPos>& info)
{
	using namespace std;
	using vector_t = typename TAAPos::ValueType;

	Grid& grid = info.grid();
	Face* nbrFaces[2];
	vector<Edge*> edges;

	TAAPos& aaPos = info.position_accessor();

	while(info.candidates_left()){
		Edge* e = info.pop_candidate();

	//	we only perform swaps on regular manifolds.
		if(GetAssociatedFaces(nbrFaces, grid, e, 2) == 2){
		//	make sure that both neighbors are triangles
			if(nbrFaces[0]->num_vertices() != 3 || nbrFaces[1]->num_vertices() != 3)
				continue;

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
					if((edges[j] != e) && !(info.is_segment(edges[j])
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
							 const vector3 &areaNormal,
							 number smallsq = SMALL_SQ)
{
	number t;

	vector3 lineDir, edgeDir, planeNormal;
	VecSubtract(lineDir, lineTo, lineFrom);
	VecNormalize(lineDir, lineDir);
	VecSubtract(edgeDir, edgeVrt2, edgeVrt1);
	number threshold = smallsq * VecLengthSq(edgeDir);
	VecNormalize(edgeDir, edgeDir);
	VecCross(planeNormal, edgeDir, areaNormal);
	if(RayPlaneIntersection(vOut, t, lineFrom, lineDir, edgeVrt1, planeNormal)){
		if(t > 0 && (VecDistanceSq(lineFrom, vOut) < VecDistanceSq(lineFrom, lineTo) + threshold)){
			return true;
		}
	}
	return false;
}



template <typename TAAPos>
bool QualityGridGeneration(Grid& grid, DelaunayInfo<TAAPos>& info,
						   number minAngle, int maxSteps)
{
	EnableDelaunayDebugSave(false);
	minAngle = max<number>(minAngle, 0);
	minAngle = min<number>(minAngle, 60);

	// maxSteps = 20;
	// UG_LOG("DEBUG: Setting maxSteps to " << maxSteps << "\n");

	using namespace std;

	using vector_t = typename TAAPos::ValueType;
	using DI = DelaunayInfo<TAAPos>;

	TAAPos aaPos = info.position_accessor();

	int stepCount = 0;

//	helper to collect neighbors
	Face* nbrFaces[2];
	queue<Vertex*> qvrts; // used during splits
	vector<Vertex*> vrts;
	vector<Edge*> edges;
	vector<Edge*> closeEdges; // used during splits
	vector<Face*> faces;

	MakeDelaunay(info);

	if(minAngle > 0 && maxSteps != 0){
		const number offCenterThresholdAngle = 1.1 * minAngle;
		const number offCenterATan = atan(deg_to_rad(0.5 * offCenterThresholdAngle));
		// UG_LOG("off-center atan: " << offCenterATan << endl);

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

		//	calculate triangle-circumcenter
			vector_t vpos[3] = {aaPos[f->vertex(0)], aaPos[f->vertex(1)], aaPos[f->vertex(2)]};
			vector_t cc;

			if(!TriangleCircumcenter(cc, vpos[0], vpos[1], vpos[2])){
				UG_LOG("Couldn't calculate triangle-circumcenter. Expect unexpected results!\n");
				UG_LOG("triangle: " << faceCenter << "\n");
				//SaveGridToFile(grid, "delaunay_debug.ugx");
				return false;
				//continue;
			}

		//	Test whether an off-center would be more appropriate ('Alper Üngörs Off-Centers')
			if(offCenterATan > 0)
			{
				number edgeLenSq[3] = {VecDistanceSq(vpos[0], vpos[1]),
									   VecDistanceSq(vpos[1], vpos[2]),
									   VecDistanceSq(vpos[2], vpos[0])};
				size_t shortestEdge = 0;
				for(size_t iedge = 1; iedge < 3; ++iedge){
					if(edgeLenSq[iedge] < edgeLenSq[shortestEdge])
						shortestEdge = iedge;
				}

				size_t vrtInd[2] = {shortestEdge, (shortestEdge + 1) % 3};
				vector_t dir[2];
				for(size_t ivrt = 0; ivrt < 2 ; ++ivrt){
					VecSubtract(dir[ivrt], vpos[vrtInd[ivrt]], cc);
				}

				number angle = rad_to_deg(VecAngle(dir[0], dir[1]));

				if(angle < offCenterThresholdAngle){
					vector_t base;
					VecScaleAdd(base, 0.5, vpos[vrtInd[0]], 0.5, vpos[vrtInd[1]]);
					vector_t ndir;
					VecSubtract(ndir, cc, base);
					VecNormalize(ndir, ndir);
					VecScale(ndir, ndir, 0.5 * sqrt(edgeLenSq[shortestEdge]) / offCenterATan);
					// UG_LOG("Adjusting cc from: " << cc << " to: ");
					VecAdd(cc, base, ndir);
					// UG_LOG(cc << endl);
				}
			}

		//	locate the triangle which contains cc. Do this by traversing edges
		//	as required. Note that since the delaunay property holds, we're only
		//	traversing edges in a circle, which does not contain any vertices.
			Edge* lastTraversedEdge = nullptr;
			Face* curFace = f;
			//vector_t startPos = CalculateCenter(f, aaPos);
			vector_t rayDir;
			VecSubtract(rayDir, cc, faceCenter);
			Vertex* pointInserted = nullptr;

			// UG_LOG("curFace: " << faceCenter << "\n");
			while(pointInserted == nullptr){
				//UG_LOG("curTri: " << CalculateCenter(curFace, aaPos) << "\n");
			//	to make things as robust as possible, we'll always intersect the
			//	same line with each edge
				Edge* nextEdge = nullptr;
				bool split = false;

				CollectAssociated(edges, grid, curFace);
				for(size_t i = 0; i < edges.size(); ++i){
					Edge* e = edges[i];
					if(e == lastTraversedEdge)
						continue;

					vector_t a;
					if(DelaunayLineLineIntersection(a, faceCenter, cc,
							aaPos[e->vertex(0)], aaPos[e->vertex(1)], faceNormal, SMALL_SQ))
					{
					//	this edge will be traversed next.
						nextEdge = e;
						number threshold =	VecDistanceSq(faceCenter, cc)
											- SMALL_SQ * EdgeLengthSq(e, aaPos);

					//	check whether e has to be split, to avoid bad aspect ratios
						split = (VecDistanceSq(faceCenter, a) > threshold);
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
					const bool isSegment = info.is_segment(nextEdge);
					split |= isSegment;
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
						vector_t center = CalculateCenter(nextEdge, aaPos);

						Vertex* vrt0 = nextEdge->vertex(0);
						Vertex* vrt1 = nextEdge->vertex(1);
						Vertex* edgeVrts[2] = {vrt0, vrt1};
						number radiusSq = VecDistanceSq(center, aaPos[vrt0]);
						number radius = sqrt(radiusSq);
						pointInserted = SplitEdge<RegularVertex>(grid, nextEdge, false);

						if(isSegment){
						//	depending on the marks of the corners of nextEdge,
						//	we may have to place the new point at a circular shell
							int dartInd = -1;
							for(int ivrt = 0; ivrt < 2; ++ivrt){
								if(info.mark(edgeVrts[ivrt]) == DI::DART){
									dartInd = ivrt;
									break;
								}
							}
							if(dartInd != -1){
								Vertex* dartVrt = edgeVrts[dartInd];
								Vertex* otherVrt = edgeVrts[(dartInd + 1) % 2];
								typename DI::Mark m = info.mark(otherVrt);
								if((m == DI::NEW_SEGMENT) || (m == DI::SHELL)){
								//	the new vertex has to be placed on a circular shell!
									info.set_mark(pointInserted, DI::SHELL);
									number dist = VecDistance(center, aaPos[dartVrt]);
									number csCur = 1;
									number csNext;

									if(dist >= 1){
										csNext = 2;
										while(csNext < dist){
											csCur = csNext;
											csNext *= 2.;
										}
									}
									else{
										csNext = 0.5;
										while(csNext > dist){
											csCur = csNext;
											csNext *= 0.5;
										}
									}

									vector_t dir;
									VecSubtract(dir, center, aaPos[dartVrt]);
									VecNormalize(dir, dir);

									if(fabs(dist - csCur) < fabs(dist - csNext))
										VecScale(dir, dir, csCur);
									else
										VecScale(dir, dir, csNext);

									VecAdd(center, aaPos[dartVrt], dir);
								
								//	we have to adjust the critical radius
									radiusSq = min(VecDistanceSq(center, aaPos[vrt0]),
												   VecDistanceSq(center, aaPos[vrt1]));
									radius = sqrt(radiusSq);
								}
							}

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
											if(info.is_segment(e) && !EdgeContains(e, pointInserted))
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

										//	if the vertrex was created during remeshing
										//	and lies in the circle, then we'll push it to vrts
											if((info.mark(vcon) == DI::NEW_INNER)
												&& (VecDistanceSq(aaPos[vcon], center) < radiusSq))
											{
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
								Face* parent = nullptr;
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
				//...
				//	But first we have to check whether an edge of the triangle
				//	would be encroached.
				//	If this is the case, we'll insert the new vertex at the triangle's
				//	center instead of the circum-center, to avoid skinny triangles
					CollectAssociated(edges, grid, curFace);
					for(size_t i = 0; i < edges.size(); ++i){
						Edge* e = edges[i];
						if(info.is_segment(e)){
							vector_t eCenter = CalculateCenter(e, aaPos);
							number eLenSq = EdgeLengthSq(e, aaPos);
							if(VecDistanceSq(eCenter, cc) < eLenSq / 4.){
								cc = CalculateCenter(curFace, aaPos);
								break;
							}
						}
					}

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
						if(!info.is_inner(faces[i_face]))
							continue;

						CollectAssociated(edges, grid, faces[i_face]);
						for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
							Edge* e = edges[i_edge];
							if(info.is_candidate(e) || info.is_segment(e))
								continue;

							Face* nbrs[2];
							if(GetAssociatedFaces(nbrs, grid, e, 2) == 2){
								if(info.is_inner(nbrs[0])
								   && info.is_inner(nbrs[1]))
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
