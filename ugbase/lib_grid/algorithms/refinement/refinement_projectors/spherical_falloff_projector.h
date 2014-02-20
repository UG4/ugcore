// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 20, 2014

#ifndef __H__UG__spherical_falloff_projector__
#define __H__UG__spherical_falloff_projector__

#include "../refinement_callbacks.h"

namespace ug{

///	\addtogroup lib_grid_algorithms_refinement
///	@{
////////////////////////////////////////////////////////////////////////
///	calculates new positions of vertices by projecting on a sphere
/**	In the range between innerRadius and outerRadius, the projection
 * gradually moves from spherical-projection to normal linear interpolation.
 *
 * Make sure to initialize the callback correctly. Use the same grid
 *	on which the refinement-operations will be performed. Make sure
 *	that aPos (given in the constructor) is attached to the vertices
 *	of the grid.
 *
 *	An uninitialized refinement-callback may not be used during refinement.
 */
template <class TAPosition>
class SphericalFalloffProjector : public IRefinementCallback
{
	public:
		SphericalFalloffProjector();

	///	make sure that aPos is attached to the vertices of the grid.
		SphericalFalloffProjector(Grid& grid, TAPosition& aPos,
								 const typename TAPosition::ValueType& center,
								 number innerRadius, number outerRadius);

		virtual ~SphericalFalloffProjector();

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
		number										m_innerRadius;
		number										m_outerRadius;
};
/// @}

}// end of namespace

#include "spherical_falloff_projector_impl.h"

#endif
