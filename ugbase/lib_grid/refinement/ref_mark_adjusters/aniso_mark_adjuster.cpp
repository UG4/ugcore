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

#include "aniso_mark_adjuster.h"

namespace ug{
	
void AnisoMarkAdjuster::
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
		const int anisoMark = ref.get_aniso_mark(f);

		if((refMark != RM_NONE) && (refMark < RM_REFINE) && anisoMark){
			grid.associated_elements_sorted(assEdges, f);
			for(size_t iedge = 0; iedge < assEdges.size(); ++iedge){
				Edge* e = assEdges[iedge];
				const RefinementMark edgeMark = ref.get_mark(e);

				if(		(anisoMark & (1<<iedge))
					&&	(edgeMark != RM_REFINE))
				{
					ref.mark(e, RM_REFINE);
				}
				
				if(e->is_constraining() && (edgeMark == RM_REFINE)){
					ref.mark(f, RM_REFINE);
					break;
				}
			}
		}
	}


	// FaceDescriptor fd;

	// for(size_t ivol = 0; ivol < vols.size(); ++ivol){
	// 	Volume* vol = vols[ivol];
	// 	const int volAnisoMark = ref.get_aniso_mark(vol);
	// 	if(!volAnisoMark || !(ref.get_mark(v) & RM_ANISOTROPIC))
	// 		continue;

	// 	const size_t numFaces = vol->num_faces();
	// 	for(size_t iface = 0; iface < numFaces; ++iface){
	// 		int sideMark = 0;
	// 		vol->face_desc(iface, fd);
	// 		const size_t numFaceEdges = fd.num_vertices();
	// 		for(size_t i = 0; i < numFaceEdges; ++i){
	// 			const int volEdgeIndex = vol->get_face_edge_index(iface, i);
	// 			sideMark |= ((volAnisoMark >> volEdgeIndex) & 1) << i;
	// 		}

	// 		Face* f = grid.get_face(&fd);
	// 		if(ref.get_mark(f) == RM_REFINE){
	// 			ref.mark(vol, RM_REFINE);
	// 			break;
	// 		}

	// 		const int oldAnisoMark = ref.get_aniso_mark(f);

	// 		if(oldAnisoMark != sideMark){
	// 			if((oldAnisoMark & sideMark) == oldAnisoMark){
	// 			//	oldAnisoMark is contained in sideMark
	// 				ref.mark_aniso(f, sideMark);
	// 			}
	// 			else{

	// 			}
	// 		}

	// 		if(		(ref.get_mark(f) == RM_REFINE)
	// 			||	(oldAnisoMark && (oldAnisoMark != sideMark)))
	// 		{
	// 			ref.mark(vol, RM_REFINE);
	// 			break;
	// 		}
	// 		else
	// 			ref.mark_aniso(f, sideMark);
	// 	}
	// }
}

}//	end of namespace
