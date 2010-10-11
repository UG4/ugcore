//	created by Martin Stepniewski, Sebastian Reiter
//	mastep@gmx.de, s.b.reiter@googlemail.com
//	y08 m12 d07

#ifndef __H__UG__SUBDIVISION_LOOP__
#define __H__UG__SUBDIVISION_LOOP__

#include "lib_grid/lg_base.h"
#include "subdivision_rules_piecewise_loop.h"

namespace ug
{

struct VrtInfo{
		VrtInfo(VertexBase* v, size_t val) : vrt(v), creaseValence(val) {}
		VertexBase* vrt;
		size_t creaseValence;//0->no crease
	};
	
////////////////////////////////////////////////////////////////////////
template <class TAVrtPos> void
ProjectToLimitPLoop(Grid& grid, TAVrtPos aPos, TAVrtPos aProjPos)
{
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
	
//todo: remove the struct and the vector. This can be performed on the fly
//		if surface-marks are supplied.
	std::vector<VrtInfo>	vrtInfos;
	
	
//	iterate through all vertices
	for(VertexBaseIterator iter = grid.vertices_begin();
		iter != grid.vertices_end(); ++iter)
	{
		VertexBase* v = *iter;
		
	//	check whether the vertex is a crease vertex or not
	//todo: this has to be more flexible
		if(IsBoundaryVertex2D(grid, v)){
		//	project the crease vertex
			EdgeBase* nbrs[2];
			size_t numNbrs = 0;
			for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(v);
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
				vector3& p0 = aaPos[GetConnectedVertex(nbrs[0], v)];
				vector3& p1 = aaPos[GetConnectedVertex(nbrs[1], v)];
				vector3 w = subdiv.proj_crease_weights();
				VecScaleAdd(aaProjPos[v], w.x, aaPos[v], w.y, p0, w.z, p1);
			}
			else
				aaProjPos[v] = aaPos[v];
		}
		else{
		//	project an inner vertex
		//	we have to calculate the valence and
		//	we have to check whether the vertex is a neighbour to a crease.
		//todo: This could be given my some kind of mark.
			bool creaseNbr = false;
			
		//	calculate the valence
		//todo: This should be performed in a common method.
			vrtInfos.clear();
			for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(v);
				iter != grid.associated_edges_end(v); ++iter)
			{
				VertexBase* nbrVrt = GetConnectedVertex(*iter, v);
			//	we have to check whether nbrVrt is a crease edge. If it is,
			//	we have to calculate its valence.
				if(IsBoundaryVertex2D(grid, v)){
					creaseNbr = true;
					size_t creaseValence = 0;
					for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(nbrVrt);
						iter != grid.associated_edges_end(nbrVrt); ++iter)
					{
						++creaseValence;
					}
					
					vrtInfos.push_back(VrtInfo(nbrVrt, creaseValence));
				}
				else
					vrtInfos.push_back(VrtInfo(nbrVrt, 0));
			}

		//	if the vertex is a crease-nbr, we'll apply special weights.
		//	if not, normal loop-projection-masks are used.
		//todo: add special weights
			{
			//	default loop projection
				VecScale(aaProjPos[v], aaPos[v],
						subdiv.proj_inner_center_weight(vrtInfos.size()));
				
				for(size_t i = 0; i < vrtInfos.size(); ++i){
					VecScaleAdd(aaProjPos[v], 1.0, aaProjPos[v],
								subdiv.proj_inner_nbr_weight(vrtInfos.size()),
								aaPos[vrtInfos[i].vrt]);
				}
			}
		}
	}
	
//	clean up
	if(usingTmpAttachment){
	//todo: swap attachment buffers
	//	detach it
		grid.detach_from_vertices(aProjPos);
	}
}

////////////////////////////////////////////////////////////////////////
//	ProjectToLimitLoop
/// projects surface vertices to their limit subdivision surface position
bool ProjectToLimitLoop(Grid& grid, APosition& aProjPos);

}//	end of namespace

#endif
