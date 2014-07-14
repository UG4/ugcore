// created by mstepnie
// martin.stepniewski@gcsc.uni-frankfurt.de
// Juli 14, 2014

#ifndef __H__UG__smooth_tetgrid_subdivision_projectors__
#define __H__UG__smooth_tetgrid_subdivision_projectors__

#include "../refinement_callbacks.h"

namespace ug{

///	Refines tetrahedral grids with Schaefer's (1) hybrid tetrahedral/octahedral scheme.
/**	Inner vertices are positioned using the Schaefer's scheme, boundary vertices
 *  are positioned using Loop's scheme.
 *
 *  Please note that the refiner works on different position attachments:
 *  A source and a target position attachment. If you use the refiner in the
 *  context of multigrid refinement, you can use the same attachment for both.
 *  However, if the callback is used to refine a flat grid, then it's crucial
 *  to use different attachments.
 *
 *  (1) Schaefer et al. 2004 - "Smooth Subdivision of Tetrahedral Meshes"
 */
template <class TAPosition>
class SubdivisionVolumesProjector : public RefinementCallbackLinear<TAPosition>
{
	private:
		typedef RefinementCallbackLinear<TAPosition> BaseClass;
		typedef typename TAPosition::ValueType		pos_type;

	public:
		using BaseClass::new_vertex;

	public:
		SubdivisionVolumesProjector();

	///	make sure that aPos and aTargetPos are attached to the vertices of the grid.
	/**	For flat grids it is crucial that aPos and aTargetPos differ.
	 *  If g is a multigrid, that aPos and aTargetPos can safely be the same.*/
		SubdivisionVolumesProjector(Grid& g,
										  TAPosition& aPos,
										  TAPosition& aTargetPos);

		virtual ~SubdivisionVolumesProjector();

		virtual void new_vertex(Vertex* vrt, Vertex* parent);
		virtual void new_vertex(Vertex* vrt, Edge* parent);
		virtual void new_vertex(Vertex* vrt, Face* parent);
		virtual void new_vertex(Vertex* vrt, Volume* parent);

	protected:
		virtual bool is_crease_vertex(Vertex* vrt);
		virtual bool is_crease_edge(Edge* edge);

	protected:
		Grid::VertexAttachmentAccessor<TAPosition>	m_aaTargetPos;
};
/// @}

}// end of namespace

#include "smooth_tetgrid_subdivision_projectors_impl.h"

#endif
