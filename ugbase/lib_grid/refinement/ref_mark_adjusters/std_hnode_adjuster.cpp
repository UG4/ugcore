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

#include "./std_hnode_adjuster.h"
#include "lib_grid/tools/periodic_boundary_manager.h"
#include "lib_grid/algorithms/debug_util.h"

namespace ug{

// marks geometric object e for refinement if it is periodic
template <typename TElem>
static void mark_if_periodic(IRefiner& ref, TElem* e) {
	if(!ref.grid())
		return;
	if(!ref.grid()->has_periodic_boundaries())
		return;

	PeriodicBoundaryManager& pbm = *ref.grid()->periodic_boundary_manager();

	// ensure element is periodic
	if(!pbm.is_periodic(e))
		return;

	RefinementMark refMark = ref.get_mark(e);
	if(pbm.is_master(e))
	{
		using SlaveContainer = typename PeriodicBoundaryManager::Group<TElem>::SlaveContainer;
		using SlaveIterator = typename PeriodicBoundaryManager::Group<TElem>::SlaveIterator;
		SlaveContainer& slaves = *pbm.slaves(e);
		for (SlaveIterator iter = slaves.begin(); iter != slaves.end(); ++iter)
			ref.mark(*iter, refMark);
	}
	else { // is slave
		ref.mark(pbm.master(e), refMark);
	}
}

void StdHNodeAdjuster::
ref_marks_changed(IRefiner& ref,
			   	  const std::vector<Vertex*>& vrts,
			   	  const std::vector<Edge*>& edges,
			   	  const std::vector<Face*>& faces,
			   	  const std::vector<Volume*>& vols)
{
	if(!ref.grid())
		return;
	Grid& grid = *ref.grid();

	Grid::edge_traits::secure_container		assEdges;
	Grid::face_traits::secure_container		assFaces;
	Grid::volume_traits::secure_container 	assVols;

////////////////////////////////
//	VERTICES
	if(node_dependency_order_1_enabled()){
	//	make sure that a hanging node is never constrained by an element which
	//	is has constrained vertices as corners
		for(size_t i_vrt = 0; i_vrt < vrts.size(); ++i_vrt){
			Vertex* vrt = vrts[i_vrt];
			if(!vrt->is_constrained())
				continue;
			ConstrainedVertex* hv = dynamic_cast<ConstrainedVertex*>(vrt);
			if(!hv)
				continue;

		//	make sure that all parents are marked
			GridObject* co = hv->get_constraining_object();
			if(co){
				if(auto e = dynamic_cast<Edge*>(co)){
					ref.mark(e);
				}
				else if(auto f = dynamic_cast<Face*>(co)){
					ref.mark(f);
				}
			}
		}
	}

////////////////////////////////
//	EDGES
	for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
		Edge* e = edges[i_edge];

		if(ref.get_mark(e) != RM_REFINE)
			continue;

	//	check whether hangingNodeOrder1 is enabled. If so, we have to check
	//	for associated hanging vertices and push them to qHVrts.
		if(node_dependency_order_1_enabled()){
			for(size_t i = 0; i < 2; ++i){
				if(e->vertex(i)->is_constrained()){
					Vertex* v = e->vertex(i);
					ref.mark(v);
				}
			}
		}

	//	depending on the type of the edge, we have to perform different operations
		if(auto cde = dynamic_cast<ConstrainedEdge*>(e))
		{
		//	the edge is a constrained edge. Make sure that its constraining edge
		//	or face will be refined.
			if(auto cge = dynamic_cast<ConstrainingEdge*>( cde->get_constraining_object()))
			{
				ref.mark(cge);
			}
			else if(auto cgf = dynamic_cast<ConstrainingFace*>( cde->get_constraining_object()))
			{
				ref.mark(cgf);
			}
			else{
				UG_ASSERT(grid.is_parallel(),
						  "Constrained edge doesn't have a constraining edge. In "
						  "a serial environment this should always be the case!"
						  << " At: " << GetGridObjectCenter(grid, e));
			}
		}
		else if(auto cge = dynamic_cast<ConstrainingEdge*>(e))
		{
		//	associated faces and volumes have to be marked
			if(grid.num_faces() > 0){
				grid.associated_elements(assFaces, cge);
				for(size_t i = 0; i < assFaces.size(); ++i){
					if (!(ref.get_mark(assFaces[i]) & RM_ANISOTROPIC)) // do not mark RM_REFINE if already marked anisotropic
						ref.mark(assFaces[i]);
				}
			}

			if(grid.num_volumes() > 0){
				grid.associated_elements(assVols, cge);
				for(size_t i = 0; i < assVols.size(); ++i){
					if (!(ref.get_mark(assVols[i]) & RM_ANISOTROPIC)) // do not mark RM_REFINE if already marked anisotropic
						ref.mark(assVols[i]);
				}
			}
		}
//NOTE: The check below was intended to replace the one above to reduce element refinement.
//		However, since constrained edges may be marked after this check has
//		been performed for a constraining edge, at least in parallel environments,
//		this implementation leads to problems (segfaults!)
		// else if(ConstrainingEdge* cge = dynamic_cast<ConstrainingEdge*>(e))
		// {
		// //	if one of the constrained objects is marked, then all associated faces
		// //	and volumes have to be marked, too
		// 	const size_t nce = cge->num_constrained_edges();
		// 	bool oneIsMarked = false;
		// 	for(size_t i = 0; i < nce; ++i){
		// 		const RefinementMark rm = ref.get_mark(cge->constrained_edge(i));
		// 		if(rm & RM_REFINE){
		// 			oneIsMarked = true;
		// 			break;
		// 		}
		// 	}

		// 	if(oneIsMarked){
		// 		if(grid.num_faces() > 0){
		// 			grid.associated_elements(assFaces, cge);
		// 			for(size_t i = 0; i < assFaces.size(); ++i){
		// 				if (!(ref.get_mark(assFaces[i]) & RM_ANISOTROPIC)) // do not mark RM_REFINE if already marked anisotropic
		// 					ref.mark(assFaces[i]);
		// 			}
		// 		}

		// 		if(grid.num_volumes() > 0){
		// 			grid.associated_elements(assVols, cge);
		// 			for(size_t i = 0; i < assVols.size(); ++i){
		// 				if (!(ref.get_mark(assVols[i]) & RM_ANISOTROPIC)) // do not mark RM_REFINE if already marked anisotropic
		// 					ref.mark(assVols[i]);
		// 			}
		// 		}
		// 	}
		// }
	}

////////////////////////////////
//	FACES
	for(size_t i_face = 0; i_face < faces.size(); ++i_face){
		Face* f = faces[i_face];
		RefinementMark refMark = ref.get_mark(f);
	//	check whether hangingNodeOrder1 is enabled. If so, we have to check
	//	for associated hanging vertices and push them to qHVrts.
		if(node_dependency_order_1_enabled()){
			for(size_t i = 0; i < f->num_vertices(); ++i){
				if(f->vertex(i)->is_constrained())
					ref.mark(f->vertex(i));
			}
		}

	//	we have to make sure that associated edges are marked.
		if(refMark != RM_LOCAL){
			grid.associated_elements(assEdges, f);
			for(size_t i = 0; i < assEdges.size(); ++i){
				if(refMark > ref.get_mark(assEdges[i]))
					ref.mark(assEdges[i], refMark);
			}
		}

	//	constrained and constraining faces require special treatment
		if(auto cdf = dynamic_cast<ConstrainedFace*>(f)){
		//	make sure that its constraining face will be refined
			if(auto cgf = dynamic_cast<ConstrainingFace*>( cdf->get_constraining_object()))
			{
				ref.mark(cgf);
			}
			else{
				UG_ASSERT(grid.is_parallel(),
						  "Constrained face doesn't have a constraining face. In "
						  "a serial environment this should always be the case!"
						  << " At: " << GetGridObjectCenter(grid, f));
			}
		}
		else if(ConstrainingFace* cgf = dynamic_cast<ConstrainingFace*>(f)){
		//	associated volumes have to be marked
			if(grid.num_volumes() > 0){
				grid.associated_elements(assVols, cgf);
				for(size_t i = 0; i < assVols.size(); ++i){
				//todo:	also check whether the local mark matches
					if (!(ref.get_mark(assVols[i]) & RM_ANISOTROPIC)) // do not mark RM_REFINE if already marked anisotropic
						ref.mark(assVols[i]);
				}
			}
		}
	}

////////////////////////////////
//	VOLUMES
	for(size_t i_vol = 0; i_vol < vols.size(); ++i_vol){
		Volume* v = vols[i_vol];
		RefinementMark refMark = ref.get_mark(v);
	//	we have to make sure that all associated edges and faces are marked.
		if(refMark != RM_LOCAL){
			grid.associated_elements(assEdges, v);
			for(size_t i = 0; i < assEdges.size(); ++i){
				if(refMark > ref.get_mark(assEdges[i]))
					ref.mark(assEdges[i], refMark);
			}

			grid.associated_elements(assFaces, v);
			for(size_t i = 0; i < assFaces.size(); ++i){
				if(refMark > ref.get_mark(assFaces[i]))
					ref.mark(assFaces[i], refMark);
			}
		}
	}

////////////////////////////////
// Periodic boundaries
	if(grid.has_periodic_boundaries()){
		for(size_t i_vrt = 0; i_vrt < vrts.size(); ++i_vrt)
			mark_if_periodic(ref, vrts[i_vrt]);

		for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge)
			mark_if_periodic(ref, edges[i_edge]);

		for(size_t i_face = 0; i_face < faces.size(); ++i_face)
			mark_if_periodic(ref, faces[i_face]);

		// omit volumes, as these are not meant to be periodic
	}
}
}// end of namespace
