/*
 * Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
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

#include "local_mark_adjuster.h"

namespace ug{
	
void LocalMarkAdjuster::
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
	
	for(size_t iface = 0; iface < faces.size(); ++iface){
		Face* f = faces[iface];
		const RefinementMark refMark = ref.get_mark(f);
		const int localMark = ref.get_local_mark(f);
		if((refMark == RM_LOCAL) && localMark){
			grid.associated_elements_sorted(assEdges, f);
			for(size_t iedge = 0; iedge < assEdges.size(); ++iedge){
				Edge* e = assEdges[iedge];
				const RefinementMark edgeMark = ref.get_mark(e);

				if(		(localMark & (1<<iedge))
					&&	(edgeMark != RM_FULL))
				{
					ref.mark(e, RM_FULL);
				}
				
				// if(e->is_constraining() && (edgeMark == RM_FULL)){
				// 	ref.mark(f, RM_FULL);
				// 	break;
				// }
			}
		}
	}


//todo: iterate over faces and check whether volume-aniso-marks and
//		associated face-aniso-marks mismatch

	Grid::face_traits::secure_container		assFaces;
	std::vector<int>	vinds;
	vinds.reserve(4);

	for(size_t ivol = 0; ivol < vols.size(); ++ivol){
		Volume* vol = vols[ivol];
		const int volLocalMark = ref.get_local_mark(vol);
		if(!volLocalMark || !(ref.get_mark(vol) & RM_LOCAL))
			continue;

		const size_t numEdges = vol->num_edges();
		for(size_t iedge = 0; iedge < numEdges; ++iedge){
			if(!(volLocalMark & (1<<iedge)))
				continue;

			Edge* e = grid.get_edge(vol, iedge);
			const RefinementMark edgeMark = ref.get_mark(e);

			if(edgeMark != RM_FULL)
				ref.mark(e, RM_FULL);
		}

		grid.associated_elements(assFaces, vol);
		for(size_t iface = 0; iface < assFaces.size(); ++iface){
			Face* f = assFaces[iface];
			Face::ConstVertexArray vrts = f->vertices();
			const size_t numVrts = f->num_vertices();

			vinds.resize(numVrts);
			for(size_t i = 0; i < numVrts; ++i){
				vinds[i] = GetVertexIndex(vol, vrts[i]);
			}


			int sideMark = 0;

			for(size_t i = 0; i < numVrts; ++i){
				const int edgeInd =
						vol->get_edge_index_from_vertices(
								vinds[i], vinds[(i+1)%numVrts]);
				
				sideMark |= ((volLocalMark >> edgeInd) & 1) << i;
			}

			if(sideMark){
				const int curMark = ref.get_mark(f);
				if(curMark == RM_FULL)
					continue;
				else if(curMark < RM_LOCAL)
					ref.mark_local(f, sideMark);

				int curLocal = ref.get_local_mark(f);

				if(curLocal != sideMark){
					if((curLocal & sideMark) == curLocal){
					//	curLocal is contained in sideMark
						ref.mark_local(f, sideMark);
					}
					else if((curLocal & sideMark) != sideMark){
					//	we have to fully refine the face, since aniso-marks do not match
						ref.mark(f, RM_FULL);
					}
				}
			}
		}
	}
}

}//	end of namespace
