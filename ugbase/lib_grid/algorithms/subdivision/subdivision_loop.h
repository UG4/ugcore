/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Sebastian Reiter, Martin Stepniewski
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

#ifndef __H__UG__SUBDIVISION_LOOP__
#define __H__UG__SUBDIVISION_LOOP__

#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "subdivision_rules_piecewise_loop.h"

namespace ug
{

///	\addtogroup lib_grid_algorithms_refinement_subdivision
///	@{

////////////////////////////////////////////////////////////////////////
///	projects all vertices in the given grid to their limit-positions using the piecewise loop scheme.
template <typename TAVrtPos> void
ProjectToLimitPLoop(Grid& grid, TAVrtPos aPos, TAVrtPos aProjPos)
{
//	position type
using pos_type = typename TAVrtPos::ValueType;
	
//	access the subdivision weights
	SubdivRules_PLoop& subdiv = SubdivRules_PLoop::inst();
	
//	if aPos and aProjPos are equal, we'll have to create a temporary
//	attachment
	bool usingTmpAttachment = false;
	if(aPos == aProjPos){
	//	create temporary attachment
		usingTmpAttachment = true;
		aProjPos = TAVrtPos();
	}
	
//	attach aProjPos if required
	if(!grid.has_vertex_attachment(aProjPos))
		grid.attach_to_vertices(aProjPos);
		
//	attachment accessors
	Grid::VertexAttachmentAccessor<TAVrtPos> aaPos(grid, aPos);
	Grid::VertexAttachmentAccessor<TAVrtPos> aaProjPos(grid, aProjPos);
	
//	vectors to hold temporary results
	using NbrInfo = SubdivRules_PLoop::NeighborInfo;
	std::vector<NbrInfo> nbrInfos;
	std::vector<Vertex*> vrts;
	std::vector<number> nbrWgts;
	
//	if volumes are contained in the grid, we currently only perform projection
//	of boundary elements (no creases...)
	if(grid.num<Volume>() > 0){
		for(VertexIterator iter = grid.vertices_begin();
			iter != grid.vertices_end(); ++iter)
		{
			Vertex* vrt = *iter;
		//	collect all surface neighbors of vrt in vrts
			vrts.clear();
			for(auto iter = grid.associated_edges_begin(vrt);
			    iter != grid.associated_edges_end(vrt); ++iter)
			{
				if(IsBoundaryEdge3D(grid, *iter))
					vrts.push_back(GetConnectedVertex(*iter, vrt));
			}
			
		//	default loop projection
			VecScale(aaProjPos[vrt], aaPos[vrt],
					subdiv.proj_inner_center_weight(vrts.size()));

			for(size_t i = 0; i < vrts.size(); ++i){
				VecScaleAdd(aaProjPos[vrt], 1.0, aaProjPos[vrt],
							subdiv.proj_inner_nbr_weight(vrts.size()),
							aaPos[vrts[i]]);
			}
		}
	}
	else{
	//	iterate through all vertices
		for(VertexIterator iter = grid.vertices_begin();
			iter != grid.vertices_end(); ++iter)
		{
			Vertex* v = *iter;
			
		//	check whether the vertex is a crease vertex or not
		//todo: this has to be more flexible
			if(IsBoundaryVertex2D(grid, v)){
			//	project the crease vertex
				Edge* nbrs[2];
				size_t numNbrs = 0;
				for(auto iter = grid.associated_edges_begin(v);
				    iter != grid.associated_edges_end(v); ++iter)
				{
					if(IsBoundaryEdge2D(grid, *iter)){
						nbrs[numNbrs] = *iter;
						++numNbrs;
						if(numNbrs == 2)
							break;
					}
				}

				if(numNbrs == 2){
					pos_type& p0 = aaPos[GetConnectedVertex(nbrs[0], v)];
					pos_type& p1 = aaPos[GetConnectedVertex(nbrs[1], v)];
					vector3 w = subdiv.proj_crease_weights();
					VecScaleAdd(aaProjPos[v], w.x(), aaPos[v], w.y(), p0, w.z(), p1);
				}
				else
					aaProjPos[v] = aaPos[v];
			}
			else{
			//	project an inner vertex
			//	we have to calculate the valence and
			//	we have to check whether the vertex is a neighbor to a crease.
			//todo: This could be given my some kind of mark.
				bool creaseNbr = false;
				
			//	collect all neighbors in an ordered set.
			//todo: the order is only important if the node is indeed a crease neighbor.
				if(!CollectSurfaceNeighborsSorted(vrts, grid, v)){
					UG_LOG("WARNING in ProjectToLimitPLoop: surface is not regular.");
					UG_LOG(" Ignoring vertex...\n");
					aaProjPos[v] = aaPos[v];
					continue;
				}
				
				nbrInfos.resize(vrts.size());
				for(size_t i = 0; i < vrts.size(); ++i)
				{
					Vertex* nbrVrt = vrts[i];
				//	we have to check whether nbrVrt is a crease edge. If it is,
				//	we have to calculate its valence.
					if(IsBoundaryVertex2D(grid, nbrVrt)){
						creaseNbr = true;
						size_t creaseValence = 0;
						for(auto iter = grid.associated_edges_begin(nbrVrt);
						    iter != grid.associated_edges_end(nbrVrt); ++iter)
						{
							++creaseValence;
						}

						nbrInfos[i] = NbrInfo(nbrVrt, creaseValence);
					}
					else
						nbrInfos[i] = NbrInfo(nbrVrt, 0);
				}

			//	if the vertex is a crease-nbr, we'll apply special weights.
			//	if not, normal loop-projection-masks are used.

				if(creaseNbr){
					number cntrWgt;
					nbrWgts.resize(vrts.size());
					subdiv.proj_inner_crease_nbr_weights(cntrWgt, &nbrWgts.front(),
														 &nbrInfos.front(), nbrInfos.size());

				//	special crease neighbor projection
					VecScale(aaProjPos[v], aaPos[v], cntrWgt);

					for(size_t i = 0; i < nbrWgts.size(); ++i){
						VecScaleAdd(aaProjPos[v], 1.0, aaProjPos[v],
									nbrWgts[i], aaPos[nbrInfos[i].nbr]);
					}
				}
				else{
				//	default loop projection
					VecScale(aaProjPos[v], aaPos[v],
							subdiv.proj_inner_center_weight(nbrInfos.size()));

					for(size_t i = 0; i < nbrInfos.size(); ++i){
						VecScaleAdd(aaProjPos[v], 1.0, aaProjPos[v],
									subdiv.proj_inner_nbr_weight(nbrInfos.size()),
									aaPos[nbrInfos[i].nbr]);
					}
				}
			}
		}
	}
	
//	clean up
	if(usingTmpAttachment){
	//	swap attachment buffers
		aaPos.swap(aaProjPos);
		
	//	detach it
		grid.detach_from_vertices(aProjPos);
	}
}

////////////////////////////////////////////////////////////////////////
///	projects all boundary vertices in the given grid to their limit-positions using the piecewise loop scheme.
template <typename TAVrtPos> void
ProjectToLimitSubdivBoundary(Grid& grid, TAVrtPos aPos, TAVrtPos aProjPos)
{
//	position type
	using pos_type = typename TAVrtPos::ValueType;
	
//	access the subdivision weights
	SubdivRules_PLoop& subdiv = SubdivRules_PLoop::inst();
	
//	if aPos and aProjPos are equal, we'll have to create a temporary
//	attachment
	bool usingTmpAttachment = false;
	if(aPos == aProjPos){
	//	create temporary attachment
		usingTmpAttachment = true;
		aProjPos = TAVrtPos();
	}
	
//	attach aProjPos if required
	if(!grid.has_vertex_attachment(aProjPos))
		grid.attach_to_vertices(aProjPos);
		
//	attachment accessors
	Grid::VertexAttachmentAccessor<TAVrtPos> aaPos(grid, aPos);
	Grid::VertexAttachmentAccessor<TAVrtPos> aaProjPos(grid, aProjPos);
	
//	iterate through all vertices
	for(VertexIterator iter = grid.vertices_begin();
		iter != grid.vertices_end(); ++iter)
	{
		Vertex* v = *iter;
		
	//	check whether the vertex is a crease vertex or not
	//todo: this has to be more flexible
		if(IsBoundaryVertex2D(grid, v)){
		//	project the crease vertex
			Edge* nbrs[2];
			size_t numNbrs = 0;
			for(auto iter = grid.associated_edges_begin(v);
			    iter != grid.associated_edges_end(v); ++iter)
			{
				if(IsBoundaryEdge2D(grid, *iter)){
					nbrs[numNbrs] = *iter;
					++numNbrs;
					if(numNbrs == 2)
						break;
				}
			}
			
			if(numNbrs == 2){
				pos_type& p0 = aaPos[GetConnectedVertex(nbrs[0], v)];
				pos_type& p1 = aaPos[GetConnectedVertex(nbrs[1], v)];
				vector3 w = subdiv.proj_crease_weights();
				VecScaleAdd(aaProjPos[v], w.x(), aaPos[v], w.y(), p0, w.z(), p1);
			}
			else
				aaProjPos[v] = aaPos[v];
		}
		else{
		//	inner vertices are not moved
			aaProjPos[v] = aaPos[v];
		}
	}
	
//	clean up
	if(usingTmpAttachment){
	//	swap attachment buffers
		aaPos.swap(aaProjPos);
		
	//	detach it
		grid.detach_from_vertices(aProjPos);
	}
}

////////////////////////////////////////////////////////////////////////
//	ProjectToLimitLoop
/// projects surface vertices to their limit subdivision surface position
bool ProjectToLimitLoop(Grid& grid, APosition& aProjPos);

/// @}	//	end of add_to_group

}//	end of namespace

#endif
