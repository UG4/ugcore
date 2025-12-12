/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#include "triangle_fill_sweep_line.h"

#include <list>
#include <vector>
#include <cassert>
#include <map>
#include <stack>
#include <algorithm>

#include "lib_grid/lg_base.h"

using namespace std;

namespace ug {

struct SweepLineEdge;

enum SweepLineVertexStatus{
	SLVS_NONE = 0,
	SLVS_START,
	SLVS_END,
	SLVS_REGULAR,
	SLVS_SPLIT,
	SLVS_MERGE
};

enum SweepLineEdgeStatus{
	SLES_UNKNOWN = 0,
	SLES_RIM
};

/**	Please note that the connections array has to be maintained manually.
 *	This is important to reduce unnecessary overhead.
 */
struct SweepLineVertex
{
	SweepLineVertex() : m_status(SweepLineVertexStatus::SLVS_NONE), isRimVertex(false)	{}
	SweepLineVertex(int ind, vector2* ptr) :
		vrtInd(ind), vrtPtr(ptr), m_status(SweepLineVertexStatus::SLVS_NONE), isRimVertex(false)	{}

	int vrtInd;
	const vector2* vrtPtr;
	unsigned char m_status;
	std::vector<SweepLineEdge*> connections;
	bool isRimVertex;
};

/**	Please note that edges do not automatically register them selves at
 *	their vertices connections array. This would make problems with
 *	temporary edge objects and copy-construction.*/
struct SweepLineEdge
{
	SweepLineEdge(SweepLineVertex* v1, SweepLineVertex* v2) :
		m_status(SweepLineEdgeStatus::SLES_UNKNOWN), m_v1(v1), m_v2(v2), m_helper(nullptr), m_numConnectedPolygons(0)
	{}

	SweepLineVertex* get_connected(const SweepLineVertex* v) const
	{
		if(m_v1 == v)
			return m_v2;
		else if(m_v2 == v)
			return m_v1;
		return nullptr;
	}

	vector2 direction() const
	{
		return vector2(m_v2->vrtPtr->x() - m_v1->vrtPtr->x(),
						m_v2->vrtPtr->y() - m_v1->vrtPtr->y());
	}

	bool contains(const SweepLineVertex* v) const {
		return (m_v1 == v) || (m_v2 == v);
	}

	unsigned char m_status;
	SweepLineVertex* m_v1;
	SweepLineVertex* m_v2;
	SweepLineVertex* m_helper;
	int m_numConnectedPolygons;
};


using SweepLineEdgeList = list<SweepLineEdge>;
using SweepLineEdgeIter = SweepLineEdgeList::iterator;
using MapEdgeCuts = multimap<number, SweepLineEdge*>;


inline bool cmp_slv(const SweepLineVertex& v1, const SweepLineVertex& v2)
{
	if(v1.vrtPtr->y() > v2.vrtPtr->y())
		return true;
	else if(v1.vrtPtr->y() == v2.vrtPtr->y())
		if(v1.vrtPtr->x() < v2.vrtPtr->x())
			return true;
	return false;
}

///	returns a value between -1 and 1, determining the magnitude of the right bend.
/**	(-1, 0): left bend,
 *	0: no bend,
 *	(0, 1): right bend
 */
number CalculateRightBendVal(vector2& dirOut,
							 number& normDotOut,
							 SweepLineEdge& e,
							 SweepLineVertex* firstVrt,
							 vector2 lastDir,
							 bool lastDirNormalized = false)
{
	if(!lastDirNormalized)
		VecNormalize(lastDir, lastDir);

//	the direction has to point away from lastVrt
	dirOut = e.direction();
	if(e.m_v1 != firstVrt)
		VecScale(dirOut, dirOut, -1.);

	VecNormalize(dirOut, dirOut);

//	automatically normalized
	vector2 lastNorm(lastDir.y(), -lastDir.x());

//	calculate the value for this configuration.
//	the result will be between -1 and 1.
//	The higher the better.
	number dDir = VecDot(lastDir, dirOut);
	normDotOut = VecDot(lastNorm, dirOut);
	number val;
	if(normDotOut <= 0)
		val = 0.5 * dDir - 0.5;
	else
		val = 0.5 - 0.5 * dDir;
	return val;
}

/**
 *	For each entry in vrtsIn an entry in vrtsOut is created.
 *	However, vertices in vrtsOut are sorted from top to bottom.
 *	Each vertex in vrtsOut references all edges in edgesOut, to
 *	which it is connected.
 *	Directions of edges in edgesOut are CCW, but only if sorting
 *	was possible (this is always the case for the outer rim).
 */
bool CreateSweepLineStructs(vector<SweepLineVertex>& vrtsOut,
							SweepLineEdgeList& edgesOut,
							const vector<vector2>& vrtsIn,
							const vector<int>& edgesIn)
{
	if(vrtsIn.size() < 3)
		return false;

	vrtsOut.resize(vrtsIn.size());

//	fill entries
	for(size_t i = 0; i < vrtsIn.size(); ++i){
		vrtsOut[i].vrtInd = i;
		vrtsOut[i].vrtPtr = &vrtsIn[i];
	}

//	sort the sweeplinevertices
	sort(vrtsOut.begin(), vrtsOut.end(), cmp_slv);

//	redirect all connection indices for the SweepLineVertices
	{
	//	an array for index indirections.
	//	associates the index of a normal vertex with the index of a slv
		vector<int> slvInd(vrtsOut.size());
		for(size_t i = 0; i < vrtsOut.size(); ++i){
			slvInd[vrtsOut[i].vrtInd] = i;
		}

	//	add connections
		for(size_t i = 0; i < edgesIn.size(); i+=2){
			edgesOut.push_back(SweepLineEdge(&vrtsOut[slvInd[edgesIn[i]]],
								&vrtsOut[slvInd[edgesIn[i+1]]]));
			vrtsOut[slvInd[edgesIn[i]]].connections.push_back(&edgesOut.back());
			vrtsOut[slvInd[edgesIn[i+1]]].connections.push_back(&edgesOut.back());
		}
	}

//	now make sure that the outer rim is in ccw order
//	start at the first vertex (this is the top vertex)
//	for each vertex follow the edge, that makes the sharpest turn to the right.
//	relative to the last checked edge. the initial direction is left. This is fine,
//	since we start at the top.
	SweepLineVertex* comingFromVrt = nullptr;
	vector2 lastDir(-1, 0);
	SweepLineVertex* lastVrt = &vrtsOut[0];
	lastVrt->m_status = SweepLineVertexStatus::SLVS_START;
	lastVrt->isRimVertex = true;

//	the start vertex has to have at least one connection
	if(lastVrt->connections.empty()){
	//	theres nothing we can do.
		return false;
	}

	while(true){
	//	find the edge that makes the sharpest turn to the right (even if it is a left turn)
	//	normally a vertex has to have two connections. However, we make an exception for the first one.
		if((lastVrt != &vrtsOut[0]) && (lastVrt->connections.size() < 2)){
		//	we're at the end of the road - this means the outer rim is open.
		//	close it and break
			edgesOut.emplace_back(lastVrt, &vrtsOut[0]);
			lastVrt->connections.push_back(&edgesOut.back());
			vrtsOut[0].connections.push_back(&edgesOut.back());
			break;
		}

		number bestVal = -2;
		vector2 nextDir(0, 0);
		SweepLineEdge* nextEdge = nullptr;
		number bestNormDot = 0;
		for(size_t i = 0; i < lastVrt->connections.size(); ++i){
			SweepLineEdge* tmpEdge= lastVrt->connections[i];
			if(tmpEdge->get_connected(lastVrt) == comingFromVrt)
				continue;
			
			vector2 tDir;
			number dNorm;
			number val = CalculateRightBendVal(tDir, dNorm, *tmpEdge, lastVrt, lastDir, true);

			if(val > bestVal){
				bestVal = val;
				nextDir = tDir;
				nextEdge = tmpEdge;
				bestNormDot  = dNorm;
			}
		}

		assert((nextEdge != nullptr) && "a new direction should have been found!");

	//	make sure that nextIter points into the right direction
		if(nextEdge->m_v1 != lastVrt)
			swap(nextEdge->m_v1, nextEdge->m_v2);

	//	based on lastDir and nextDir we can assign the new status to lastVrt, if none has been
	//	previously assigned.
		if(lastVrt->m_status == SweepLineVertexStatus::SLVS_NONE){
			if(lastDir.y() > 0){
				if(nextDir.y() < 0){
					if(bestNormDot < 0)
						lastVrt->m_status = SweepLineVertexStatus::SLVS_START;
					else
						lastVrt->m_status = SweepLineVertexStatus::SLVS_SPLIT;
				}
				else if(nextDir.y() > 0)
					lastVrt->m_status = SweepLineVertexStatus::SLVS_REGULAR;
				else{
					if(nextDir.x() < 0)
						lastVrt->m_status = SweepLineVertexStatus::SLVS_REGULAR;
					else
						lastVrt->m_status = SweepLineVertexStatus::SLVS_SPLIT;
				}

			//	if an inner edge leaves lastVrt upwards, then we have to
			//	regard it as a regular vertex.
				if(lastVrt->m_status == SweepLineVertexStatus::SLVS_SPLIT){
					for(size_t i = 0; i < lastVrt->connections.size(); ++i){
						if(cmp_slv(*lastVrt->connections[i]->get_connected(lastVrt), *lastVrt)){
							lastVrt->m_status = SweepLineVertexStatus::SLVS_REGULAR;
							break;
						}
					}
				}
			}
			else if(lastDir.y() < 0){	//	lastDir.y() < 0 from now on
				if(nextDir.y() > 0){
					if(bestNormDot > 0)
						lastVrt->m_status = SweepLineVertexStatus::SLVS_MERGE;
					else
						lastVrt->m_status = SweepLineVertexStatus::SLVS_END;
				}
				else if(nextDir.y() < 0)
					lastVrt->m_status = SweepLineVertexStatus::SLVS_REGULAR;
				else{
					if(nextDir.x() > 0)
						lastVrt->m_status = SweepLineVertexStatus::SLVS_REGULAR;
					else
						lastVrt->m_status = SweepLineVertexStatus::SLVS_MERGE;
				}

			//	if an inner edge leaves lastVrt downwards, then we have to
			//	regard it as a regular vertex.
				if(lastVrt->m_status == SweepLineVertexStatus::SLVS_MERGE){
					for(size_t i = 0; i < lastVrt->connections.size(); ++i){
						if(!cmp_slv(*lastVrt->connections[i]->get_connected(lastVrt), *lastVrt)){
							lastVrt->m_status = SweepLineVertexStatus::SLVS_REGULAR;
							break;
						}
					}
				}
			}
			else{ // lastDir.y() == 0 from now on
				if(lastDir.x() < 0){
					if(nextDir.y() >= 0)
						lastVrt->m_status = SweepLineVertexStatus::SLVS_REGULAR;
					else
						lastVrt->m_status = SweepLineVertexStatus::SLVS_START;
				}
				else{
					if(nextDir.y() <= 0)
						lastVrt->m_status = SweepLineVertexStatus::SLVS_REGULAR;
					else
						lastVrt->m_status = SweepLineVertexStatus::SLVS_END;
				}
			}
		}

	//	mark the current edge as a rim edge
		nextEdge->m_status = SweepLineEdgeStatus::SLES_RIM;

	//	get the next vertex
		comingFromVrt = lastVrt;
		lastVrt = nextEdge->m_v2;
		lastVrt->isRimVertex = true;
		lastDir = nextDir;
	//	if the vertex has already been processed, we have to quit here
		if(lastVrt->m_status != SweepLineVertexStatus::SLVS_NONE)
			break;
	}

//	now make sure that vertices that do not lie on the outer rim have a correct status
	for(size_t i = 0; i < vrtsOut.size(); ++i){
		SweepLineVertex& vrt = vrtsOut[i];
		if(vrt.m_status == SweepLineVertexStatus::SLVS_NONE){
		//	check the position of connected vertices
			bool gotUpper = false;
			bool gotLower = false;
			bool gotEqualYLeft = false;
			bool gotEqualYRight = false;
			bool gotRight = false;
			for(size_t j = 0; j < vrt.connections.size(); ++j){
				SweepLineVertex& connVrt = *(vrt.connections[j]->get_connected(&vrt));

				if(connVrt.vrtPtr->x() > vrt.vrtPtr->x())
					gotRight = true;

				if(connVrt.vrtPtr->y() < vrt.vrtPtr->y())
					gotLower = true;
				else if(connVrt.vrtPtr->y() > vrt.vrtPtr->y())
					gotUpper = true;
				else{
					if(connVrt.vrtPtr->x() > vrt.vrtPtr->x())
						gotEqualYRight = true;
					else
						gotEqualYLeft = true;
				}
			}

		//	if there is only one connection, we have to treat the vertex specially
			if(vrt.connections.size() == 1){
				if(gotUpper)
					vrt.m_status = SweepLineVertexStatus::SLVS_MERGE;
				else if(gotLower)
					vrt.m_status = SweepLineVertexStatus::SLVS_SPLIT;
				else{
					if(gotRight)
						vrt.m_status = SweepLineVertexStatus::SLVS_SPLIT;
					else
						vrt.m_status = SweepLineVertexStatus::SLVS_MERGE;
				}
			}
			else{ // more than one connection from now on
				if(gotUpper && gotLower)
					vrt.m_status = SweepLineVertexStatus::SLVS_REGULAR;
				else{
					if(gotEqualYLeft && gotEqualYRight)
						vrt.m_status = SweepLineVertexStatus::SLVS_REGULAR;
					else if(gotUpper){
						if(gotEqualYRight)
							vrt.m_status = SweepLineVertexStatus::SLVS_REGULAR;
						else
							vrt.m_status = SweepLineVertexStatus::SLVS_MERGE;
					}
					else if(gotLower){
						if(gotEqualYLeft)
							vrt.m_status = SweepLineVertexStatus::SLVS_REGULAR;
						else
							vrt.m_status = SweepLineVertexStatus::SLVS_SPLIT;
					}
					else
						vrt.m_status = SweepLineVertexStatus::SLVS_REGULAR;
				}
			}

			assert((vrt.m_status != SLVS_NONE) && "no status found");
			if(vrt.m_status == SweepLineVertexStatus::SLVS_NONE){
				UG_LOG("ERROR in CreateSweeplineStructs: Could not assign status.\n");
				return false;
			}
		}
	}

//	we're finally through with it!!!
	return true;
}

void PrintDebugInfos(const vector<SweepLineVertex>& vrts,
					 const list<SweepLineEdge>& edges)
{
	UG_LOG("SweepLine debug infos:");
	UG_LOG("  num vrts: " << vrts.size() << endl);
	UG_LOG("  vertex types (top to bottom): ");
	for(size_t i = 0; i < vrts.size(); ++i){
		if(vrts[i].isRimVertex){
			UG_LOG("rim-");
		}
		else{
			UG_LOG("inner-");
		}

		switch(vrts[i].m_status){
			case SweepLineVertexStatus::SLVS_NONE:		UG_LOG("none"); break;
			case SweepLineVertexStatus::SLVS_START:	UG_LOG("start"); break;
			case SweepLineVertexStatus::SLVS_END:		UG_LOG("end"); break;
			case SweepLineVertexStatus::SLVS_REGULAR:	UG_LOG("regular"); break;
			case SweepLineVertexStatus::SLVS_SPLIT:	UG_LOG("split"); break;
			case SweepLineVertexStatus::SLVS_MERGE:	UG_LOG("merge"); break;
			default:			UG_LOG("unknown"); break;
		}
		UG_LOG(", ");
	}

//	log inner edges
	UG_LOG("\n  inner edges:\n");
	for(list<SweepLineEdge>::const_iterator iter = edges.begin(); iter != edges.end(); ++iter){
		if(iter->m_status != SweepLineEdgeStatus::SLES_RIM)
		{
			UG_LOG("     - " << *iter->m_v1->vrtPtr << ", " << *iter->m_v2->vrtPtr << endl);
		}
	}
	UG_LOG(endl);

	UG_LOG("  connections:\n");
	for(size_t i = 0; i < vrts.size(); ++i){
		const SweepLineVertex& v = vrts[i];
		UG_LOG("    vertex " << *v.vrtPtr << " is connected to: ");
		for(size_t j = 0; j < v.connections.size(); ++j){
			UG_LOG(*v.connections[j]->get_connected(&v)->vrtPtr << ", ");
		}
		UG_LOG("\n");
	}

}

bool SweepLineEdgeIntersectsSweepLine(number& xOut,
									  const SweepLineEdge& edge,
									  number sweepLineY)
{
//	first check whether both points lie above or below the sweepline
	if((edge.m_v1->vrtPtr->y() < sweepLineY)
	   && (edge.m_v2->vrtPtr->y() < sweepLineY))
	   return false;

	if((edge.m_v1->vrtPtr->y() >= sweepLineY)
	   && (edge.m_v2->vrtPtr->y() >= sweepLineY))
	   return false;

	vector2 dir = edge.direction();
	if(dir.y() != 0){
		number s = (sweepLineY - edge.m_v1->vrtPtr->y()) / dir.y();
		//if((s >= 0) && (s <= 1.)){
		if(s < 0) s = 0;
		if(s > 1.) s = 1.;

		xOut = edge.m_v1->vrtPtr->x() + s * dir.x();
		return true;
	}

	return false;
}

SweepLineEdge* GetEdgeOnTheLeft(MapEdgeCuts& edgeCuts, SweepLineVertex& v)
{
//	find the edge in edgeCuts which is directly left of v
	SweepLineEdge* leftEdge = nullptr;
	number leftEdgeVal = 0;

	for(auto iter = edgeCuts.begin(); iter != edgeCuts.end(); ++iter)
	{
		if(iter->first < v.vrtPtr->x()){
			if(!iter->second->contains(&v)){
				bool takeNew = true;
			//	if a second edge with identical cut exists (this can happen
			//	if two edges are cut at their upper vertex), we have to make
			//	sure to return the edge, which points farther to the right.
				if(leftEdge){
					if(leftEdgeVal == iter->first){
					//	check whether iter->second points farther to the right.
						vector2 dirOld = leftEdge->direction();
						vector2 dirNew = iter->second->direction();
						if(dirOld.y() > 0)
							VecScale(dirOld, dirOld, -1);
						if(dirNew.y() > 0)
							VecScale(dirNew, dirNew, -1);

						VecNormalize(dirOld, dirOld);
						VecNormalize(dirNew, dirNew);

						if(dirOld.x() > dirNew.x())
							takeNew = false;
					}
				}

				if(takeNew){
					leftEdge = iter->second;
					leftEdgeVal = iter->first;
				}
			}
		}
		else
			break;
	}
	return leftEdge;
}

bool EdgeExists(SweepLineVertex* v0, SweepLineVertex* v1)
{
	for(size_t i = 0; i < v0->connections.size(); ++i){
		if(v0->connections[i]->contains(v1)){
			return true;
		}
	}
	return false;
}

///	inserts new edges so that edges contains a set of monotone polynomes
bool SweepLine_CreateMonotones(vector<SweepLineVertex>& vrts,
								list<SweepLineEdge>& edges)
{
//	create sets of monotone polynomes.
	MapEdgeCuts edgeCuts;

	for(size_t i_vrt = 0; i_vrt < vrts.size(); ++i_vrt){
		SweepLineVertex& v = vrts[i_vrt];

		number sweepLineY = v.vrtPtr->y();
	//	update the edgeCutsMap.
		MapEdgeCuts	tmap;
		for(auto iter = edgeCuts.begin(); iter != edgeCuts.end(); ++iter)
		{
		//	if the helper was set to nullptr, we'll ignore the edge
			if(iter->second->m_helper){
			//	make sure that the edge still intersects the sweepLine.
				number x;
				if(SweepLineEdgeIntersectsSweepLine(x, *iter->second, sweepLineY))
					tmap.insert(make_pair(x, iter->second));
			}
		}
//only for debugging!
//vector2 tmpPos = *v.vrtPtr;

		edgeCuts.swap(tmap);

		switch(v.m_status){
			case SweepLineVertexStatus::SLVS_NONE:
				break;
			case SweepLineVertexStatus::SLVS_START:
			//	rim-edges that leave v and that have v as their first
			//	vertex, have to be added to edgeCuts.

			//	add all edges connected to the start vertex to edgeCuts.
			//	Note that at least one edge with has exterior to its right
			//	will be added.
				for(size_t i = 0; i < v.connections.size(); ++i){
					edgeCuts.insert(make_pair(v.vrtPtr->x(), v.connections[i]));
					v.connections[i]->m_helper = &v;
					/*
					if((v.connections[i]->m_status == SLES_RIM)
					   && (v.connections[i]->m_v1 == &v))
					{
						edgeCuts.insert(make_pair(v.vrtPtr->x(), v.connections[i]));
						v.connections[i]->m_helper = &v;
					}
					else if(v.connections[i]->m_status != SLES_RIM){
					//	inner edges that go down have to be added, too
						if(v.connections[i]->get_connected(&v)->vrtPtr->y() < v.vrtPtr->y()){
							edgeCuts.insert(make_pair(v.vrtPtr->x(), v.connections[i]));
							v.connections[i]->m_helper = &v;
						}
					}*/
				}
				break;
			case SweepLineVertexStatus::SLVS_END:
				{
				//	get the incoming edge
					SweepLineEdge* incoming = nullptr;
					for(size_t i = 0; i < v.connections.size(); ++i){
						if((v.connections[i]->m_status == SweepLineEdgeStatus::SLES_RIM)
						   && (v.connections[i]->m_v2 == &v))
						{
							incoming = v.connections[i];
							break;
						}
					}

					//assert((incoming != nullptr) && "An incoming edge has to exist!");
					if(!incoming){
						UG_LOG("ERROR in SweepLine_CreateMonotones: couldn't find incoming-edge.\n");
						return false;
					}

					if(incoming->m_helper != nullptr){
						if(incoming->m_helper->m_status == SweepLineVertexStatus::SLVS_MERGE){
						//	insert a diagonal between helper and v
							edges.emplace_back(incoming->m_helper, &v);
							incoming->m_helper->connections.push_back(&edges.back());
							v.connections.push_back(&edges.back());
						}
					}

				//	remove incoming from the set of edge-cuts.
				//	This can be done by setting its helper to nullptr.
				//	The rest will be handled later on
					incoming->m_helper = nullptr;
				}break;
			case SweepLineVertexStatus::SLVS_REGULAR:
				{
				//	get the direction of the incoming edge
					if(v.isRimVertex && (v.connections.size() == 2)){
						SweepLineEdge* incoming = nullptr;
						for(size_t i = 0; i < v.connections.size(); ++i){
							if((v.connections[i]->m_status == SweepLineEdgeStatus::SLES_RIM)
							   && (v.connections[i]->m_v2 == &v))
							{
								incoming = v.connections[i];
								break;
							}
						}

						//assert((incoming != nullptr) && "An incoming edge has to exist!");

						if(!incoming){
							UG_LOG("ERROR in SweepLine_CreateMonotones: couldn't find incoming-edge.\n");
							return false;
						}

						vector2 dir = incoming->direction();
						bool areaIsToTheRight = false;
						if(dir.y() < 0)
							areaIsToTheRight = true;
						else if(dir.y() == 0){
							if(dir.x() > 0)
								areaIsToTheRight = true;
						}

						if(areaIsToTheRight){
							if(!incoming->m_helper){
								UG_LOG("ERROR in SweepLine_CreateMonotones (SLVS_REGULAR-rim): incoming-edge has no helper at vertex " << vrts[i_vrt].vrtInd << ".\n");
								return false;
							}
							if(incoming->m_helper->m_status == SweepLineVertexStatus::SLVS_MERGE){
							//	insert a diagonal between helper and v
								edges.emplace_back(incoming->m_helper, &v);
								incoming->m_helper->connections.push_back(&edges.back());
								v.connections.push_back(&edges.back());

							//	delete incoming from edgeCuts
								incoming->m_helper = nullptr;
							}

						//	add the outgoing edge to edgeCuts
							SweepLineEdge* outgoing = nullptr;
							for(size_t i = 0; i < v.connections.size(); ++i){
								if((v.connections[i]->m_status == SweepLineEdgeStatus::SLES_RIM)
								   && (v.connections[i]->m_v1 == &v))
								{
									outgoing = v.connections[i];
									break;
								}
							}

							//assert((outgoing != nullptr) && "An outgoing edge has to exist!");
							if(!outgoing){
								UG_LOG("ERROR in SweepLine_CreateMonotones: couldn't find outgoing-edge.\n");
								return false;
							}

							edgeCuts.insert(make_pair(v.vrtPtr->x(), outgoing));
							outgoing->m_helper = &v;
						}
						else{ //	area is on the left
						//	find the edge in edgeCuts which is directly left of v
							SweepLineEdge* leftEdge = GetEdgeOnTheLeft(edgeCuts, v);

							if(leftEdge){
								if(!leftEdge->m_helper){
									UG_LOG("ERROR in SweepLine_CreateMonotones (SLVS_REGULAR-rim): left-edge has no helper at vertex " << vrts[i_vrt].vrtInd << ".\n");
									return false;
								}
								if(leftEdge->m_helper->m_status == SweepLineVertexStatus::SLVS_MERGE){
								//	add a diagonal and replace the helper
									edges.emplace_back(leftEdge->m_helper, &v);
									leftEdge->m_helper->connections.push_back(&edges.back());
									v.connections.push_back(&edges.back());
								}
								leftEdge->m_helper = &v;
							}
						}
					}
					else{
					//	the vertex is an inner regular vertex or a rim-vertex which is connected to
					//	inner edges.
					//	check all incoming edges (edges that come from above)
					//	make sure to not check the new ones
						size_t initialSize = v.connections.size();
						for(size_t i = 0; i < initialSize; ++i){
							SweepLineEdge& incoming = *v.connections[i];
							SweepLineVertex* conn = incoming.get_connected(&v);

						//	rim vertices have to be treated specially
							if(v.isRimVertex){
							//	if the edge is a rim edge, then it has to end at v
								if((incoming.m_status == SweepLineEdgeStatus::SLES_RIM)
								   && (incoming.m_v1 == &v))
								{
									continue;
								}
							}

						//	make sure that the edge is really an incoming edge
							if(conn->vrtPtr->y() < v.vrtPtr->y())
								continue;
							else if(conn->vrtPtr->y() == v.vrtPtr->y()){
								if(conn->vrtPtr->x() > v.vrtPtr->x())
									continue;
							}

						//	ok. it is.
							//assert((incoming.m_helper != nullptr) && "a helper has to exist");
							if(incoming.m_helper)
							{
								if(incoming.m_helper->m_status == SweepLineVertexStatus::SLVS_MERGE){
								//	insert a diagonal between helper and v
									edges.emplace_back(incoming.m_helper, &v);
									incoming.m_helper->connections.push_back(&edges.back());
									v.connections.push_back(&edges.back());

								//	delete incoming from edgeCuts
									incoming.m_helper = nullptr;
								}
							}
						}

					//	now find the edge in edgeCuts, which is directly left of v
						SweepLineEdge* leftEdge = GetEdgeOnTheLeft(edgeCuts, v);

						if(leftEdge){
							if(!leftEdge->m_helper){
								UG_LOG("ERROR in SweepLine_CreateMonotones (SLVS_REGULAR): left-edge has no helper at vertex " << i_vrt << ".\n");
								return false;
							}
							if(leftEdge->m_helper->m_status == SweepLineVertexStatus::SLVS_MERGE){
							//	add a diagonal and replace the helper
								edges.emplace_back(leftEdge->m_helper, &v);
								leftEdge->m_helper->connections.push_back(&edges.back());
								v.connections.push_back(&edges.back());
							}
							leftEdge->m_helper = &v;
						}

					//	finally add all outgoing edges to edgeCuts
						for(size_t i = 0; i < v.connections.size(); ++i){
							SweepLineEdge& outgoing = *v.connections[i];
							SweepLineVertex* conn = outgoing.get_connected(&v);

						//	make sure that the edge is really an outgoing edge
							if(conn->vrtPtr->y() > v.vrtPtr->y())
								continue;
							else if(conn->vrtPtr->y() == v.vrtPtr->y()){
								if(conn->vrtPtr->x() < v.vrtPtr->x())
									continue;
							}

						//	ok. it is outgoing
							edgeCuts.insert(make_pair(v.vrtPtr->x(), v.connections[i]));
							outgoing.m_helper = &v;
						}
					}
				}break;
			case SweepLineVertexStatus::SLVS_SPLIT:
				{
				//	find the edge in edgeCuts which is directly left of v
					SweepLineEdge* leftEdge = GetEdgeOnTheLeft(edgeCuts, v);

					//assert((leftEdge != nullptr) && "a edge on the left has to exist.");

					if(!leftEdge){
						UG_LOG("ERROR in SweepLine_CreateMonotones: couldn't find left-edge.\n");
						return false;
					}

					if(!leftEdge->m_helper){
						UG_LOG("ERROR in SweepLine_CreateMonotones (SLVS_SPLIT): left-edge has no helper at vertex " << i_vrt << ".\n");
						return false;
					}

				//	create a new diagonal from leftEdges helper to v
					edges.emplace_back(leftEdge->m_helper, &v);
					leftEdge->m_helper->connections.push_back(&edges.back());
					v.connections.push_back(&edges.back());

				//	replace the helper by v
					leftEdge->m_helper = &v;
				//	add the outgoing edges of v to the map
				//	if the vertex is an inner vertex, we have to add all connected edges.
					for(size_t i = 0; i < v.connections.size(); ++i){
						if((!v.isRimVertex) || (v.connections[i]->m_v1 == &v))
						{
							edgeCuts.insert(make_pair(v.vrtPtr->x(), v.connections[i]));
							v.connections[i]->m_helper = &v;
						}
					}
				}break;
			case SweepLineVertexStatus::SLVS_MERGE:
				{
				//	we differentiate between rim vertices and inner vertices, since this gives a speedup.
					if(v.isRimVertex){
					//	get the incoming edge
						SweepLineEdge* incoming = nullptr;
						for(size_t i = 0; i < v.connections.size(); ++i){
							if((v.connections[i]->m_status == SweepLineEdgeStatus::SLES_RIM)
							   && (v.connections[i]->m_v2 == &v))
							{
								incoming = v.connections[i];
								break;
							}
						}
						//assert((incoming != nullptr) && "An incoming edge has to exist!");
						//assert((incoming->m_helper != nullptr) && "The edge has to have a helper!");

						if(!incoming){
							UG_LOG("ERROR in SweepLine_CreateMonotones: couldn't find incoming-edge.\n");
							return false;
						}

						if(!incoming->m_helper){
							UG_LOG("ERROR in SweepLine_CreateMonotones (SLVS_MERGE): incoming-edge has no helper at vertex " << i_vrt << ".\n");
							return false;
						}

						if(incoming->m_helper->m_status == SweepLineVertexStatus::SLVS_MERGE){
						//	insert a diagonal between helper and v
							edges.emplace_back(incoming->m_helper, &v);
							incoming->m_helper->connections.push_back(&edges.back());
							v.connections.push_back(&edges.back());
						}

					//	remove incoming from the set of edge-cuts.
					//	This can be done by setting its helper to nullptr.
					//	The rest will be handled later on
						incoming->m_helper = nullptr;
					}
					else{
					//	check all incoming edges
						for(size_t i = 0; i < v.connections.size(); ++i){
							SweepLineEdge* incoming = v.connections[i];
							if(incoming->m_helper){
								if(incoming->m_helper->m_status == SweepLineVertexStatus::SLVS_MERGE){
									edges.emplace_back(incoming->m_helper, &v);
									incoming->m_helper->connections.push_back(&edges.back());
									v.connections.push_back(&edges.back());
								}
							//	remove incoming from the set of edge-cuts.
							//	This can be done by setting its helper to nullptr.
							//	The rest will be handled later on
								incoming->m_helper = nullptr;
							}
						}
					}

				//	find the edge in edgeCuts which is directly left of v
					SweepLineEdge* leftEdge = GetEdgeOnTheLeft(edgeCuts, v);

					//assert((leftEdge != nullptr) && "an edge on the left has to exist.");
					//assert(leftEdge->m_helper && "edge has no helper");

					if(!leftEdge){
						UG_LOG("ERROR in SweepLine_CreateMonotones: couldn't find left-edge.\n");
						return false;
					}

					if(!leftEdge->m_helper){
						UG_LOG("ERROR in SweepLine_CreateMonotones (SLVS_MERGE): left-edge has no helper at vertex " << i_vrt << ".\n");
						return false;
					}

					if(leftEdge->m_helper != &v
					   && leftEdge->m_helper->m_status == SweepLineVertexStatus::SLVS_MERGE)
					{
					//	add a diagonal and replace the helper
						edges.emplace_back(leftEdge->m_helper, &v);
						leftEdge->m_helper->connections.push_back(&edges.back());
						v.connections.push_back(&edges.back());
					}

					leftEdge->m_helper = &v;

				}break;
		}
	}

	return true;
}

bool TriangleFill_SweepLine(vector<int>& facesOut,
							const vector<vector2>& srcVrtsOrig,
							/*const */vector<int>& srcEdges)
{
//	first check, that all edges are fine
	size_t numVrts = srcVrtsOrig.size();
	for(size_t i = 0; i < srcEdges.size(); ++i){
		if(srcEdges[i] < 0 || srcEdges[i] >= (int)numVrts){
			UG_LOG("Bad edge index " << srcEdges[i] << " at index "
					<< i << " in srcEdges in TriangleFill_SweepLine. "
					"Please check your edge-input array!\n");
			return false;
		}
	}

//	THIS IS A NASTY HACK: By casting all the values to float, we can avoid
//	some rounding issues... This has to be improved!
//	At the same time we're transforming the pointset, to improve the algorithm
//	for axis aligned geometries (which occur quite often...)
	vector<vector2> srcVrts(srcVrtsOrig.size());
	for(size_t i = 0; i < srcVrtsOrig.size(); ++i){
		srcVrts[i].x() = (float)(srcVrtsOrig[i].x());
		srcVrts[i].y() = (float)(srcVrtsOrig[i].y() - 0.333213214 * srcVrtsOrig[i].x());
	}

//	create an array with SweepLineVertices
//	be careful not to change this array after it has been set up,
//	since otherwise pointers might be invalidated.
	vector<SweepLineVertex> vrts;
	list<SweepLineEdge> edges;

	if(!CreateSweepLineStructs(vrts, edges, srcVrts, srcEdges)){
		UG_LOG("CreateSweepLineStructs failed.\n");
		//UG_LOG("Make sure not to select vertices that don't have connected edges.\n");
		UG_LOG("Make sure to select a closed outer edge-chain!\n");
		UG_LOG("Aborting.\n");
		return false;
	}

//PrintDebugInfos(vrts, edges);

	if(!SweepLine_CreateMonotones(vrts, edges)){
		UG_LOG("SweepLine_CreateMonotones failed.\n");
		UG_LOG("Make sure to select a closed outer edge-chain!\n");
		//PrintDebugInfos(vrts, edges);
		return false;
	}

/*
//DEBUG BEGIN
PrintDebugInfos(vrts, edges);
//	write all edges to srcEdges for debug purposes

srcEdges.clear();
for(SweepLineEdgeIter iter = edges.begin(); iter != edges.end(); ++iter){
	srcEdges.push_back(iter->m_v1->vrtInd);
	srcEdges.push_back(iter->m_v2->vrtInd);
}
//return false;
//DEBUG END
*/

//	triangulate monotone polygons
//	start at the top
//	while there is a candidate
//	- find its two rightmost edges that are not yet complete.
//	- if there are no, pick the next candidate and restart the iteration (or quit if none is left)
//	- there are two. Follow them until they reach a common point by always taking sharp left corners for the
//	left branch and sharp right corners for the right branch (if viewed in direction of travel.
//	edges can only be used if they are not yet complete.)
//	both branches are then stored in a monotone-polygon and the polygonCounter of each edge is increased by 1.
	int branchInd = 0;	// used for error-logging...
	for(size_t i_main = 0; i_main < vrts.size();)
	{
		SweepLineVertex& monotoneTop = vrts[i_main];
	//	find its two rightmost edges which have not yet been processed
		number bestVal[2] = {-2, -2};
		SweepLineEdge* branch[2] = {nullptr, nullptr};
		vector2 downDir(0, -1);//compare to straight down ray.
		for(size_t i = 0; i < monotoneTop.connections.size(); ++i){
			SweepLineEdge& e = *monotoneTop.connections[i];
		//	e is only a candidate for a branch, if it meets the following reqirements:
		//	- e is a rim edge and hasn't yet got a connected polygon.
		//	- e is an inner edge and has at most one connected polygon
			if(((e.m_status == SweepLineEdgeStatus::SLES_RIM) && e.m_numConnectedPolygons == 0)
			   || ((e.m_status != SweepLineEdgeStatus::SLES_RIM) && (e.m_numConnectedPolygons < 2)))
			{
				vector2 tDir;
				number dNorm;
				number val = CalculateRightBendVal(tDir, dNorm, e, &monotoneTop, downDir, true);
				if(val > bestVal[0]){
				//	move content of branch[0] to branch[1]
					bestVal[1] = bestVal[0];
					branch[1] = branch[0];
				//	set new values
					bestVal[0] = val;
					branch[0] = &e;
				}
				else if(val > bestVal[1]){
				//	set new values
					bestVal[1] = val;
					branch[1] = &e;
				}
			}
		}

	//	if we didn't find two valid branches, we have to continue with the next vertex
	//	from top to bottom. (Checking branch[1] is enough here).
		if(!branch[1]){
			++i_main;
			continue;
		}

	//	we found two. We have to increase their connectedPolygons member
		branch[0]->m_numConnectedPolygons++;
		branch[1]->m_numConnectedPolygons++;

	//	branch[0] now contains the first edge of the left branch and branch[1]
	//	the first edge of the right branch.
	//	since the branches are monotone, we can follow them from top to bottom.
	//	on the fly we increase the connectedPolygons counter of each edge.
	//	we use a stack for all encountered but unprocessed vertices
		stack<SweepLineVertex*> stk;
		stk.push(&monotoneTop);

		int lastBranchInd = -1; //0: left, 1:right
		while(true){
		//	check whether we walk on the left or on the right branch
		//	first get the lower vertex of each branch
			SweepLineVertex* vLeft = nullptr;
			SweepLineVertex* vRight = nullptr;

			if(cmp_slv(*branch[0]->m_v1, *branch[0]->m_v2))
				vLeft = branch[0]->m_v2;
			else
				vLeft = branch[0]->m_v1;

			if(cmp_slv(*branch[1]->m_v1, *branch[1]->m_v2))
				vRight = branch[1]->m_v2;
			else
				vRight = branch[1]->m_v1;

		//	the higher vertex is the next one to take
			SweepLineVertex* cur;
			int curBranchInd;
			if(cmp_slv(*vLeft, *vRight)){
				cur = vLeft;
				curBranchInd = 0;
			}
			else{
				cur = vRight;
				curBranchInd = 1;
			}

		//	if cur is equal to the top of the stack, we reached the bottom and
		//	are done. This normally shouldn't happen, since we exit at another
		//	check when searching for the next branch normally.
			//assert((cur != stk.top()) && "the loop should have already been terminated!");
			if(cur == stk.top()){
				UG_LOG("exiting since current vertex was already on the stack: " << *cur->vrtPtr);
				UG_LOG(" with index: " << cur->vrtInd << "\n");
				UG_LOG("branchIndex: " << branchInd << "\n");
				UG_LOG("ERROR in TriangleFill_SweepLine: Algorithm failed."
						" Make sure that your geometry does not contain multiple vertices"
						" at the same position (invoke RemoveDoubles on the whole selection)!\n");
				return false;
				break;
			}

		//	if there was no last branch, we'll simply insert the vertex into the
		//	stack and do some administrative stuff
			if(lastBranchInd == -1){
				stk.push(cur);
				//lastBranchInd = curBranchInd;  // never used
			}
			else{
			//	if the last branch was on the other side we can now build triangles
			//	for all stack-members
				if(lastBranchInd != curBranchInd){
				//	pop all elements from the stack and build triangles
				//	we need the first
					SweepLineVertex* v1 = stk.top(); stk.pop();
					assert((!stk.empty()) && "at least two vertices should lie in the stack at this position.");

					SweepLineVertex* oldTop = v1;//	we have to readd this later on
					while(!stk.empty()){
					//	get the next vertex
						SweepLineVertex* v2 = stk.top(); stk.pop();

					//	create a triangle with correct orientation
						if(curBranchInd == 0){
							facesOut.push_back(cur->vrtInd);
							facesOut.push_back(v1->vrtInd);
							facesOut.push_back(v2->vrtInd);
						}
						else{
							facesOut.push_back(cur->vrtInd);
							facesOut.push_back(v2->vrtInd);
							facesOut.push_back(v1->vrtInd);
						}

					//	prepare the next iteration
						v1 = v2;
					}

				//	readd the last two vertices
					stk.push(oldTop);
					stk.push(cur);
				}
				else{
				//	pop vertices from the stack and try to build triangles.
				//	if a triangle can't be build, we're done.
					SweepLineVertex* v1 = stk.top(); stk.pop();
					assert((!stk.empty()) && "at least two vertices should lie in the stack at this position.");

					while(!stk.empty()){
					//	get the next vertex
						SweepLineVertex* v2 = stk.top(); stk.pop();
					//	check whether the triangle can be build
						vector2 dir;
						VecSubtract(dir, *v1->vrtPtr, *cur->vrtPtr);
						vector2 norm(dir.y(), -dir.x());
						VecNormalize(norm, norm);
						VecSubtract(dir, *v2->vrtPtr, *v1->vrtPtr);
						VecNormalize(dir, dir);
						number dot = VecDot(dir, norm);

					//	the result has to be interpreted depending on the branch
						if(curBranchInd == 0){
							if(dot > 0.0001){
							//	create the triangle
								facesOut.push_back(cur->vrtInd);
								facesOut.push_back(v2->vrtInd);
								facesOut.push_back(v1->vrtInd);
							//	v2 is the new v1
								v1 = v2;
							}
							else{
							//	can't build triangle. restore and exit
								stk.push(v2);
								break;
							}
						}
						else{
							if(dot < -0.0001){
							//	create the triangle
								facesOut.push_back(cur->vrtInd);
								facesOut.push_back(v1->vrtInd);
								facesOut.push_back(v2->vrtInd);
							//	v2 is the new v1
								v1 = v2;
							}
							else{
							//	can't build triangle. restore and exit
								stk.push(v2);
								break;
							}
						}
					}

				//	readd v1 and cur to the stack
					stk.push(v1);
					stk.push(cur);
				}
			}

		//	now get the next branch edge
		//	if cur was on the left branch, we have to get the next edge which makes the
		//	hardest turn to the left (in direction of travel). The opposite is the case
		//	if we are on the right branch.
		//	in any case - if the next edge points upwards, we're done.
			number multiplyer = 1.0;	// left-right helper
			if(curBranchInd == 0)
				multiplyer = -1.0;

			SweepLineVertex* conn = branch[curBranchInd]->get_connected(cur);
			vector2 lastDir;
			VecSubtract(lastDir, *cur->vrtPtr, *conn->vrtPtr);
			VecNormalize(lastDir, lastDir);
			number bestVal = -2.0;
			SweepLineEdge* nextEdge = nullptr;
			vector2 newDir(0, 1);
			for(size_t i = 0; i < cur->connections.size(); ++i){
				SweepLineEdge& e = *cur->connections[i];
			//	make sure not to readd the current branch-edge.
				if(&e == branch[curBranchInd])
					continue;

				vector2 tDir;
				number dNorm;
				number val = multiplyer * CalculateRightBendVal(tDir, dNorm, e, cur, lastDir, true);

				if(val > bestVal){
				//	set new values
					bestVal = val;
					nextEdge = &e;
					newDir = tDir;
				}
			}

		//	if no new edge has been found we have to quit.
			if(!nextEdge){
				break;
			}

		//	if a new edge has been found which points upwards, we're done, too
			if(newDir.y() > 0){
				break;
			}
			else if((newDir.y() == 0) && (newDir.x() < 0)){
				break;
			}

		//	everythings fine. replace the branch edge with the new edge
			branch[curBranchInd] = nextEdge;
		//	increase the connectedPolygonCount
			branch[curBranchInd]->m_numConnectedPolygons++;
			lastBranchInd = curBranchInd;
		}

		//if(branchInd == 3)	break;
		++branchInd;
	}

	return true;
}

bool TriangleFill_SweepLine(std::vector<int>& facesOut,
							const std::vector<vector3>& srcVrts,
							/*const */std::vector<int>& srcEdges)
{
//	transform the pointset to 2d
	std::vector<vector2> vrts2d(srcVrts.size());
	TransformPointSetTo2D(&vrts2d.front(), &srcVrts.front(),
						  srcVrts.size());

	return TriangleFill_SweepLine(facesOut, vrts2d, srcEdges);
}

}//	namespace ug
