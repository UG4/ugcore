#ifndef __H__UG__loop_subdivision_projectors__
#define __H__UG__loop_subdivision_projectors__

#include "../refinement_callbacks.h"

namespace ug{

///	\addtogroup lib_grid_algorithms_refinement
///	@{

///	Applies smooth subdivision rules on boundary elements.
/**	Inner vertices are positioned using the standard linear scheme,
 *  boundary vertices are positioned using b-spline weights.
 *
 *  Please note that the refiner works on different position attachments:
 *  A source and a target position attachment. If you use the refiner in the
 *  context of multigrid refinement, you can use the same attachment for both.
 *  However, if the callback is used to refine a flat grid, then it's crucial
 *  to use different attachments..
 */
template <class TAPosition>
class SubdivisionLoopBoundaryProjector : public RefinementCallbackLinear<TAPosition>
{
	private:
		typedef RefinementCallbackLinear<TAPosition> BaseClass;
		typedef typename TAPosition::ValueType		pos_type;

	public:
		using BaseClass::new_vertex;

	public:
		SubdivisionLoopBoundaryProjector();

	///	make sure that aPos and aTargetPos are attached to the vertices of the grid.
	/**	For flat grids it is crucial that aPos and aTargetPos differ.
	 *  If g is a multigrid, that aPos and aTargetPos can safely be the same.*/
		SubdivisionLoopBoundaryProjector(Grid& g,
										 TAPosition& aPos,
										 TAPosition& aTargetPos);

		virtual ~SubdivisionLoopBoundaryProjector();

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


///	Refines grids with an extended loop scheme.
/**	Inner vertices are positioned using the loop scheme, boundary vertices
 *  are positioned using b-spline weights.
 *
 *  Please note that the refiner works on different position attachments:
 *  A source and a target position attachment. If you use the refiner in the
 *  context of multigrid refinement, you can use the same attachment for both.
 *  However, if the callback is used to refine a flat grid, then it's crucial
 *  to use different attachments.
 */
template <class TAPosition>
class SubdivisionLoopProjector : public RefinementCallbackLinear<TAPosition>
{
	private:
		typedef RefinementCallbackLinear<TAPosition> BaseClass;
		typedef typename TAPosition::ValueType		pos_type;

	public:
		using BaseClass::new_vertex;

	public:
		SubdivisionLoopProjector();

	///	make sure that aPos and aTargetPos are attached to the vertices of the grid.
	/**	For flat grids it is crucial that aPos and aTargetPos differ.
	 *  If g is a multigrid, that aPos and aTargetPos can safely be the same.*/
		SubdivisionLoopProjector(Grid& g,
										  TAPosition& aPos,
										  TAPosition& aTargetPos);

		virtual ~SubdivisionLoopProjector();

		virtual void new_vertex(Vertex* vrt, Vertex* parent);
		virtual void new_vertex(Vertex* vrt, Edge* parent);
		virtual void new_vertex(Vertex* vrt, Face* parent);
		virtual void new_vertex(Vertex* vrt, Volume* parent);

		virtual void consider_as_crease_edge(Grid::edge_traits::callback cbIsCrease);

	protected:
		virtual bool is_crease_vertex(Vertex* vrt);
		virtual bool is_crease_edge(Edge* edge);

	protected:
		Grid::VertexAttachmentAccessor<TAPosition>	m_aaTargetPos;
		Grid::edge_traits::callback					m_cbIsCrease;
};
/// @}

}// end of namespace

#include "loop_subdivision_projectors_impl.h"

#endif
