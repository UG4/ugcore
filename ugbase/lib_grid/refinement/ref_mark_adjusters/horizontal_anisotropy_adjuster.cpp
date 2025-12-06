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

#include "common/math/ugmath_types.h"
#include "horizontal_anisotropy_adjuster.h"
//ø #include "common/util/vec_for_each.h"

namespace ug{

static inline bool IsVertical(const vector1&, const vector1&)
{
	return false;
}

static inline bool IsVertical(const vector2& from, const vector2& to)
{
	return fabs(from[0] - to[0]) < SMALL;
}

static inline bool IsVertical(const vector3& from, const vector3& to)
{
	return (fabs(from[0] - to[0]) < SMALL) && (fabs(from[1] - to[1]) < SMALL);
}


static inline bool NotVertical(const vector1&, const vector1&)
{
	return true;
}

static inline bool NotVertical(const vector2& from, const vector2& to)
{
	return fabs(from[0] - to[0]) > SMALL;
}

static inline bool NotVertical(const vector3& from, const vector3& to)
{
	return (fabs(from[0] - to[0]) > SMALL) || (fabs(from[1] - to[1]) > SMALL);
}


template <typename TAPos>
void HorizontalAnisotropyAdjuster<TAPos>::
ref_marks_changed(IRefiner& ref,
				  const std::vector<Vertex*>& vrts,
				  const std::vector<Edge*>& edges,
				  const std::vector<Face*>& faces,
				  const std::vector<Volume*>& vols)
{
	UG_COND_THROW(!ref.grid(), "A grid has to be associated with the given"
					" refiner in order to adjust refinement marks.");
	
	Grid& grid = *ref.grid();
	position_accessor_t aaPos = Grid::VertexAttachmentAccessor<position_attachment_t>(grid, m_aPos);

	Grid::edge_traits::secure_container		assEdges;
	Grid::face_traits::secure_container		assFaces;
	Grid::volume_traits::secure_container	assVols;

	for(size_t _vfeI = 0; _vfeI < vols.size(); ++_vfeI){ Volume* vol = vols[_vfeI];{
		if(!(ref.get_mark(vol) & RM_ANISOTROPIC))
			continue;	// only process anisotropic volumes

		grid.associated_elements(assFaces, vol);
		for(size_t _vfeI = 0; _vfeI < assFaces.size(); ++_vfeI){ Face* f = assFaces[_vfeI];{
			if(!(ref.get_mark(f) & RM_ANISOTROPIC))
				ref.mark(f, RM_ANISOTROPIC);
		}};
	}};


	for(size_t _vfeI = 0; _vfeI < faces.size(); ++_vfeI){ Face* f = faces[_vfeI];{
		if(!(ref.get_mark(f) & RM_ANISOTROPIC))
			continue;

		grid.associated_elements(assEdges, f);
		bool horizontalMarked = false;
		bool verticalMarked = false;
		bool hasVerticalEdge = false;
		for(size_t _vfeI = 0; _vfeI < assEdges.size(); ++_vfeI){ Edge* e = assEdges[_vfeI];{
			const bool isVertical = IsVertical(aaPos[e->vertex(0)], aaPos[e->vertex(1)]);
			hasVerticalEdge |= isVertical;

			if(ref.get_mark(e) & (RM_ANISOTROPIC | RM_REFINE | RM_CLOSURE)){
				if(isVertical)
					verticalMarked = true;
				else
					horizontalMarked = true;
			}
		}};

	//	if the face contains a vertical edge and has a marked
	//	horizontal edge, we'll make sure that all horizontal
	//	edges are marked.
	//	If however a vertical edge is also marked, there's nothing we can do
		if(hasVerticalEdge && horizontalMarked && !verticalMarked){
			for(size_t _vfeI = 0; _vfeI < assEdges.size(); ++_vfeI){ Edge* e = assEdges[_vfeI];{
				if(NotVertical(aaPos[e->vertex(0)], aaPos[e->vertex(1)]))
					ref.mark(e);
			}};
		}
	}};


//	mark all associated unmarked faces and volumes of marked edges as closure elements
	if(grid.num_volumes() > 0){
		for(size_t _vfeI = 0; _vfeI < edges.size(); ++_vfeI){ Edge*e = edges[_vfeI];{
			grid.associated_elements(assVols, e);
			for(size_t _vfeI = 0; _vfeI < assVols.size(); ++_vfeI){ Volume* vol = assVols[_vfeI];{
				if(!(ref.get_mark(vol) & (RM_REFINE | RM_ANISOTROPIC | RM_CLOSURE)))	
					ref.mark(vol, RM_CLOSURE);
			}};
		}};
	}

	if(grid.num_faces() > 0){
		for(size_t _vfeI = 0; _vfeI < edges.size(); ++_vfeI){ Edge*e = edges[_vfeI];{
			grid.associated_elements(assFaces, e);
			for(size_t _vfeI = 0; _vfeI < assFaces.size(); ++_vfeI){ Face* f = assFaces[_vfeI];{
				if(!(ref.get_mark(f) & (RM_REFINE | RM_ANISOTROPIC | RM_CLOSURE)))	
					ref.mark(f, RM_CLOSURE);
			}};
		}};
	}
}

template class HorizontalAnisotropyAdjuster<AVector1>;
template class HorizontalAnisotropyAdjuster<AVector2>;
template class HorizontalAnisotropyAdjuster<AVector3>;

}//	end of namespace
