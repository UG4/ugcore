#ifndef __H__UG__cylindrical_falloff_projector__
#define __H__UG__cylindrical_falloff_projector__

#include "../refinement_callbacks.h"

namespace ug{

///	\addtogroup lib_grid_algorithms_refinement
///	@{

////////////////////////////////////////////////////////////////////////
///	calculates new positions of vertices by projecting on a cylinder
/**	Only vertices inside innerRadius are projected to a cylinder.
 * The ones outside of outerRadius are positioned through normal linear interpolation.
 * The ones in between are gradually processed from cylindrical-projection to
 * linear interpolation.
 *
 * Make sure to initialize the callback correctly. Use the same grid
 *	on which the refinement-operations will be performed. Make sure
 *	that aPos (given in the constructor) is attached to the vertices
 *	of the grid.
 *
 *	An uninitialized refinement-callback may not be used during refinement.
 */
template <class TAPosition>
class CylindricalFalloffProjector : public IRefinementCallback
{
	public:
		CylindricalFalloffProjector();

	///	make sure that aPos is attached to the vertices of the grid.
		CylindricalFalloffProjector(Grid& grid, TAPosition& aPos,
								   const typename TAPosition::ValueType& center,
								   const typename TAPosition::ValueType& axis,
								   number innerRadius, number outerRadius);

		virtual ~CylindricalFalloffProjector();

		virtual void new_vertex(Vertex* vrt, Vertex* parent);
		virtual void new_vertex(Vertex* vrt, Edge* parent);
		virtual void new_vertex(Vertex* vrt, Face* parent);
		virtual void new_vertex(Vertex* vrt, Volume* parent);

		virtual int current_pos(number* coordsOut, Vertex* vrt, int maxCoords);

	protected:
		template <class TElem>
		void perform_projection(Vertex* vrt, TElem* parent);

	protected:
		typedef typename TAPosition::ValueType		pos_type;

		Grid* 										m_pGrid;
		Grid::VertexAttachmentAccessor<TAPosition>	m_aaPos;
		pos_type									m_center;
		pos_type									m_axis;
		number										m_innerRadius;
		number										m_outerRadius;
};

/// @}

}// end of namespace

#include "cylindrical_falloff_projector_impl.h"

#endif
