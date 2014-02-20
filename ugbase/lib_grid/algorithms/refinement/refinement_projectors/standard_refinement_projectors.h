// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 20, 2014

#ifndef __H__UG__standard_refinement_projectors__
#define __H__UG__standard_refinement_projectors__

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
class RefinementCallbackSphere : public IRefinementCallback
{
	public:
		RefinementCallbackSphere();

	///	make sure that aPos is attached to the vertices of the grid.
		RefinementCallbackSphere(Grid& grid, TAPosition& aPos,
								 const typename TAPosition::ValueType& center,
								 number radius);

		virtual ~RefinementCallbackSphere();

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
		number										m_radius;
};

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
class RefinementCallbackSphericalFalloff : public IRefinementCallback
{
	public:
		RefinementCallbackSphericalFalloff();

	///	make sure that aPos is attached to the vertices of the grid.
		RefinementCallbackSphericalFalloff(Grid& grid, TAPosition& aPos,
								 const typename TAPosition::ValueType& center,
								 number innerRadius, number outerRadius);

		virtual ~RefinementCallbackSphericalFalloff();

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

////////////////////////////////////////////////////////////////////////
///	calculates new positions of vertices by projecting on a cylinder
/**	Make sure to initialize the callback correctly. Use the same grid
 *	on which the refinement-operations will be performed. Make sure
 *	that aPos (given in the constructor) is attached to the vertices
 *	of the grid.
 *
 *	An uninitialized refinement-callback may not be used during refinement.
 */
template <class TAPosition>
class RefinementCallbackCylinder : public IRefinementCallback
{
	public:
		RefinementCallbackCylinder();

	///	make sure that aPos is attached to the vertices of the grid.
		RefinementCallbackCylinder(Grid& grid, TAPosition& aPos,
								   const typename TAPosition::ValueType& center,
								   const typename TAPosition::ValueType& axis,
								   number radius);

		virtual ~RefinementCallbackCylinder();

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
		number										m_radius;
};


class RefinementCallback_IntersectCylinder : public RefinementCallbackLinear<APosition>
{
	public:
		using RefinementCallbackLinear<APosition>::new_vertex;

	public:
		RefinementCallback_IntersectCylinder();

	///	make sure that aPos is attached to the vertices of the grid.
		RefinementCallback_IntersectCylinder(Grid& grid, const vector3& center,
											 const vector3& axis, number radius,
											 APosition& aPos = aPosition);

		virtual ~RefinementCallback_IntersectCylinder();

		virtual void new_vertex(Vertex* vrt, Edge* parent);

	protected:
		vector3 m_center;
		vector3 m_axis;
		number m_radius;
};

////////////////////////////////////////////////////////////////////////
///	calculates new positions by cutting parent edges with a plane
/**	For each edge the intersection of the edge with the initially
 *	given plane is calculated and used as new point. Vertices
 *	created on other geometric objects are treated as in the linear case.
 *
 *	Make sure to initialise the callback correctly. Use the same grid
 *	on which the refinement-operations will be performed. Make sure
 *	that aPos (given in the constructor) is attached to the vertices
 *	of the grid.
 *
 *	An uninitialized refinement-callback may not be used during refinement.
 */
class RefinementCallbackEdgePlaneCut : public RefinementCallbackLinear<APosition>
{
	public:
		using RefinementCallbackLinear<APosition>::new_vertex;

	public:
		RefinementCallbackEdgePlaneCut();

	///	make sure that aPos is attached to the vertices of the grid.
		RefinementCallbackEdgePlaneCut(Grid& grid, const vector3& p,
										const vector3& n,
										APosition& aPos = aPosition);

		virtual ~RefinementCallbackEdgePlaneCut();

		virtual void new_vertex(Vertex* vrt, Edge* parent);

	protected:
		vector3 m_p;
		vector3 m_n;
};


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
class RefinementCallbackSubdivBoundary : public RefinementCallbackLinear<TAPosition>
{
	private:
		typedef RefinementCallbackLinear<TAPosition> BaseClass;
		typedef typename TAPosition::ValueType		pos_type;

	public:
		using BaseClass::new_vertex;

	public:
		RefinementCallbackSubdivBoundary();

	///	make sure that aPos and aTargetPos are attached to the vertices of the grid.
	/**	For flat grids it is crucial that aPos and aTargetPos differ.
	 *  If g is a multigrid, that aPos and aTargetPos can safely be the same.*/
		RefinementCallbackSubdivBoundary(Grid& g,
										 TAPosition& aPos,
										 TAPosition& aTargetPos);

		virtual ~RefinementCallbackSubdivBoundary();

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
class RefinementCallbackSubdivisionLoop : public RefinementCallbackLinear<TAPosition>
{
	private:
		typedef RefinementCallbackLinear<TAPosition> BaseClass;
		typedef typename TAPosition::ValueType		pos_type;

	public:
		using BaseClass::new_vertex;

	public:
		RefinementCallbackSubdivisionLoop();

	///	make sure that aPos and aTargetPos are attached to the vertices of the grid.
	/**	For flat grids it is crucial that aPos and aTargetPos differ.
	 *  If g is a multigrid, that aPos and aTargetPos can safely be the same.*/
		RefinementCallbackSubdivisionLoop(Grid& g,
										  TAPosition& aPos,
										  TAPosition& aTargetPos);

		virtual ~RefinementCallbackSubdivisionLoop();

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


#include "standard_refinement_projectors_impl.hpp"

#endif
