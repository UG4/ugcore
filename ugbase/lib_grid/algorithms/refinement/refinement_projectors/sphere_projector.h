#ifndef __H__UG__sphere_projector__
#define __H__UG__sphere_projector__

#include "../refinement_callbacks.h"

namespace ug{

///	\addtogroup lib_grid_algorithms_refinement
///	@{
////////////////////////////////////////////////////////////////////////
///	calculates new positions of vertices by projecting on a sphere
/**	Make sure to initialize the callback correctly. Use the same grid
 *	on which the refinement-operations will be performed. Make sure
 *	that aPos (given in the constructor) is attached to the vertices
 *	of the grid.
 *
 *	An uninitialized refinement-callback may not be used during refinement.
 */
template <class TAPosition>
class SphereProjector : public IRefinementCallback
{
	public:
		SphereProjector();

	///	make sure that aPos is attached to the vertices of the grid.
		SphereProjector(Grid& grid, TAPosition& aPos,
								 const typename TAPosition::ValueType& center);

		virtual ~SphereProjector();

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
};
/// @}

}// end of namespace

#include "sphere_projector_impl.h"

#endif
