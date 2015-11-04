// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 20, 2014

#ifndef __H__UG__misc_refinement_projectors__
#define __H__UG__misc_refinement_projectors__

#include "../refinement_callbacks.h"

namespace ug{

///	\addtogroup lib_grid_algorithms_refinement
///	@{

///	If a refined edge intersects the given cylinder, the new vertex will be placed at the intersection.
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




/// @}

}// end of namespace

#endif
