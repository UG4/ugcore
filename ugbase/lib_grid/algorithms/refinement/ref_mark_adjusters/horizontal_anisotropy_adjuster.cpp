// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include "common/math/ugmath_types.h"
#include "horizontal_anisotropy_adjuster.h"
#include "common/util/vec_for_each.h"

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


template <class TAPos>
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

	for_each_in_vec(Volume* vol, vols){
		if(!(ref.get_mark(vol) & RM_ANISOTROPIC))
			continue;	// only process anisotropic volumes

		grid.associated_elements(assFaces, vol);
		for_each_in_vec(Face* f, assFaces){
			if(!(ref.get_mark(f) & RM_ANISOTROPIC))
				ref.mark(f, RM_ANISOTROPIC);
		}end_for;
	}end_for;


	for_each_in_vec(Face* f, faces){
		if(!(ref.get_mark(f) & RM_ANISOTROPIC))
			continue;

		grid.associated_elements(assEdges, f);
		bool horizontalMarked = false;
		bool verticalMarked = false;
		bool hasVerticalEdge = false;
		for_each_in_vec(Edge* e, assEdges){
			const bool isVertical = IsVertical(aaPos[e->vertex(0)], aaPos[e->vertex(1)]);
			hasVerticalEdge |= isVertical;

			if(ref.get_mark(e) & (RM_ANISOTROPIC | RM_REFINE | RM_CLOSURE)){
				if(isVertical)
					verticalMarked = true;
				else
					horizontalMarked = true;
			}
		}end_for;

	//	if the face contains a vertical edge and has a marked
	//	horizontal edge, we'll make sure that all horizontal
	//	edges are marked.
	//	If however a vertical edge is also marked, there's nothing we can do
		if(hasVerticalEdge && horizontalMarked && !verticalMarked){
			for_each_in_vec(Edge* e, assEdges){
				if(NotVertical(aaPos[e->vertex(0)], aaPos[e->vertex(1)]))
					ref.mark(e);
			}end_for;
		}
	}end_for;


//	mark all associated unmarked faces and volumes of marked edges as closure elements
	if(grid.num_volumes() > 0){
		for_each_in_vec(Edge*e, edges){
			grid.associated_elements(assVols, e);
			for_each_in_vec(Volume* vol, assVols){
				if(!(ref.get_mark(vol) & (RM_REFINE | RM_ANISOTROPIC | RM_CLOSURE)))	
					ref.mark(vol, RM_CLOSURE);
			}end_for;
		}end_for;
	}

	if(grid.num_faces() > 0){
		for_each_in_vec(Edge*e, edges){
			grid.associated_elements(assFaces, e);
			for_each_in_vec(Face* f, assFaces){
				if(!(ref.get_mark(f) & (RM_REFINE | RM_ANISOTROPIC | RM_CLOSURE)))	
					ref.mark(f, RM_CLOSURE);
			}end_for;
		}end_for;
	}
}

template class HorizontalAnisotropyAdjuster<AVector1>;
template class HorizontalAnisotropyAdjuster<AVector2>;
template class HorizontalAnisotropyAdjuster<AVector3>;

}//	end of namespace
