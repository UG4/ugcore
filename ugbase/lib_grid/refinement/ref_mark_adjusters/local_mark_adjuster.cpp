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
		if(ref.marked_local(f)){
			grid.associated_elements_sorted(assEdges, f);
			const size_t numEdges = f->num_edges();
			for(size_t iedge = 0; iedge < numEdges; ++iedge){
				Edge* e = assEdges[iedge];
				if(ref.get_local_edge_mark(f, e)){
					if(!ref.marked_full(e))
						ref.mark(e, RM_FULL);
				}
			}

			if(f->is_constraining()){
			//	check which edges are also constraining and have a constrained
			//	vertex. If the local mark does not match those constraints,
			//	we have to perform a full refine.
			//todo: instead of a full refine, only regularize associated vol-marks?
				int constraintMark = 0;
				for(size_t iedge = 0; iedge < numEdges; ++iedge){
					if(ConstrainingEdge* cge =
							dynamic_cast<ConstrainingEdge*>(assEdges[iedge]))
					{
						if(cge->num_constrained_vertices())
							constraintMark |= (1<<iedge);
					}
				}

				if((constraintMark & ref.get_local_mark(f)) != constraintMark){
					ref.mark(f, RM_FULL);
				}
			}
		}

	}



	Grid::face_traits::secure_container		assFaces;
	std::vector<int>	vinds;
	vinds.reserve(4);

	for(size_t ivol = 0; ivol < vols.size(); ++ivol){
		Volume* vol = vols[ivol];
		if(!ref.marked_local(vol))
			continue;

		const size_t numEdges = vol->num_edges();
		for(size_t iedge = 0; iedge < numEdges; ++iedge){
			Edge* e = grid.get_edge(vol, iedge);
			if(!ref.get_local_edge_mark(vol, e))
				continue;

			const RefinementMark edgeMark = ref.get_mark(e);

			if(edgeMark != RM_FULL)
				ref.mark(e, RM_FULL);
		}

		grid.associated_elements(assFaces, vol);
		for(size_t iface = 0; iface < assFaces.size(); ++iface){
			Face* f = assFaces[iface];
			const int sideMark = ref.get_local_face_mark(vol, f);
			if(sideMark){
				const int curMark = ref.get_mark(f);
				if(curMark == RM_FULL)
					continue;
				else if(curMark < RM_LOCAL){
					ref.mark_local(f, sideMark);
				}
				else{
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
}

}//	end of namespace
