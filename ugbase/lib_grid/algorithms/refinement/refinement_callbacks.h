// created by Sebastian Reiter
// s.b.reiter@googlemail.com
//	y10 m07 d08

#ifndef __H__LIB_GRID__REFINEMENT_CALLBACKS__
#define __H__LIB_GRID__REFINEMENT_CALLBACKS__

#include <vector>
#include "lib_grid/lg_base.h"

namespace ug
{

///	\addtogroup lib_grid_algorithms_refinement
///	@{

////////////////////////////////////////////////////////////////////////
///	can be used to calculate new positions of vertices created during refinement.
class IRefinementCallback
{
	public:
		virtual ~IRefinementCallback()	{}
	///	called when a new vertex was created from an old vertex.
		virtual void new_vertex(VertexBase* vrt, VertexBase* parent) = 0;
	///	called when a new vertex was created from an old edge.
		virtual void new_vertex(VertexBase* vrt, EdgeBase* parent) = 0;
	///	called when a new vertex was created from an old face.
		virtual void new_vertex(VertexBase* vrt, Face* parent) = 0;
	///	called when a new vertex was created from an old volume.
		virtual void new_vertex(VertexBase* vrt, Volume* parent) = 0;

	///	callback for vertices in flat grids.
	/**	called for old vertices in flat grids which are used during refinement,
	 *	since no new vertices will be created here.
	 *	Calls new_vertex(vrt, vrt) by default.
	 *
	 *	Please note that this method won't be called during multigrid refinement.*/
		virtual void flat_grid_vertex_encountered(VertexBase* vrt)	{new_vertex(vrt, vrt);}

	///	Gives access to the position of the current vertex
	/**	The method returns the number of coordinates associated with the vertex
	 * and writes those coordinates to coordsOut. No more than maxCoords will
	 * be written to coordsOut. coordsOut thus has to be at least as big as maxCoords.*/
		virtual int current_pos(number* coordsOut, VertexBase* vrt, int maxCoords) = 0;

	protected:
	///	A helper implementation of current_pos available for derived classes.
		template <class TAttachmentAccessor>
		int current_pos_helper(number* coordsOut, VertexBase* vrt, int maxCoords,
							   TAttachmentAccessor& aaPos);
};



////////////////////////////////////////////////////////////////////////
///	calculates new positions of vertices through linear interpolation.
/**	Make sure to initialize the callback correctly. Use the same grid
 *	on which the refinement-operations will be performed. Make sure
 *	that aPos (given in the constructor) is attached to the vertices
 *	of the grid.
 *
 *	An uninitialized refinement-callback may not be used during refinement.
 */
template <class TAPosition>
class RefinementCallbackLinear : public IRefinementCallback
{
	public:
		RefinementCallbackLinear();
		
	///	make sure that aPos is attached to the vertices of the grid.
		RefinementCallbackLinear(Grid& grid, TAPosition& aPos);
	
		virtual ~RefinementCallbackLinear();
		
		virtual void new_vertex(VertexBase* vrt, VertexBase* parent);
		virtual void new_vertex(VertexBase* vrt, EdgeBase* parent);
		virtual void new_vertex(VertexBase* vrt, Face* parent);
		virtual void new_vertex(VertexBase* vrt, Volume* parent);
		
		virtual int current_pos(number* coordsOut, VertexBase* vrt, int maxCoords);

	protected:
		typedef typename TAPosition::ValueType		pos_type;
		
		Grid* 										m_pGrid;
		Grid::VertexAttachmentAccessor<TAPosition>	m_aaPos;
};


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

		virtual void new_vertex(VertexBase* vrt, VertexBase* parent);
		virtual void new_vertex(VertexBase* vrt, EdgeBase* parent);
		virtual void new_vertex(VertexBase* vrt, Face* parent);
		virtual void new_vertex(VertexBase* vrt, Volume* parent);

		virtual int current_pos(number* coordsOut, VertexBase* vrt, int maxCoords);

	protected:
		template <class TElem>
		void perform_projection(VertexBase* vrt, TElem* parent);

	protected:
		typedef typename TAPosition::ValueType		pos_type;

		Grid* 										m_pGrid;
		Grid::VertexAttachmentAccessor<TAPosition>	m_aaPos;
		pos_type									m_center;
		number										m_radius;
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

		virtual void new_vertex(VertexBase* vrt, VertexBase* parent);
		virtual void new_vertex(VertexBase* vrt, EdgeBase* parent);
		virtual void new_vertex(VertexBase* vrt, Face* parent);
		virtual void new_vertex(VertexBase* vrt, Volume* parent);

		virtual int current_pos(number* coordsOut, VertexBase* vrt, int maxCoords);

	protected:
		template <class TElem>
		void perform_projection(VertexBase* vrt, TElem* parent);

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

		virtual void new_vertex(VertexBase* vrt, EdgeBase* parent);

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
		
		virtual void new_vertex(VertexBase* vrt, EdgeBase* parent);
		
	protected:
		vector3 m_p;
		vector3 m_n;
};


////////////////////////////////////////////////////////////////////////
///	calculates new positions by offsetting new vertices to along its parents normal.
/**
 *	Make sure to initialise the callback correctly. Use the same grid
 *	on which the refinement-operations will be performed. Make sure
 *	that aPos (given in the constructor) is attached to the vertices
 *	of the grid.
 *
 *	An uninitialized refinement-callback may not be used during refinement.
 */
class RefinementCallbackFractal : public RefinementCallbackLinear<APosition>
{
	public:
		using RefinementCallbackLinear<APosition>::new_vertex;
		
	public:
		RefinementCallbackFractal();
		
	///	make sure that aPos is attached to the vertices of the grid.
		RefinementCallbackFractal(Grid& grid, number scaleFac,
								  APosition& aPos = aPosition);
	
		virtual ~RefinementCallbackFractal();
		
		virtual void new_vertex(VertexBase* vrt, EdgeBase* parent);
		virtual void new_vertex(VertexBase* vrt, Face* parent);
		
		inline void set_scale_fac(number scaleFac)	{m_scaleFac = scaleFac;}
		inline number get_scale_fac()				{return m_scaleFac;}
		
	protected:
		number m_scaleFac;
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
		
		virtual void new_vertex(VertexBase* vrt, VertexBase* parent);
		virtual void new_vertex(VertexBase* vrt, EdgeBase* parent);
		virtual void new_vertex(VertexBase* vrt, Face* parent);
		virtual void new_vertex(VertexBase* vrt, Volume* parent);
		
	protected:
		virtual bool is_crease_vertex(VertexBase* vrt);
		virtual bool is_crease_edge(EdgeBase* edge);

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
		
		virtual void new_vertex(VertexBase* vrt, VertexBase* parent);
		virtual void new_vertex(VertexBase* vrt, EdgeBase* parent);
		virtual void new_vertex(VertexBase* vrt, Face* parent);
		virtual void new_vertex(VertexBase* vrt, Volume* parent);
		
	protected:
		virtual bool is_crease_vertex(VertexBase* vrt);
		virtual bool is_crease_edge(EdgeBase* edge);

	protected:
		Grid::VertexAttachmentAccessor<TAPosition>	m_aaTargetPos;
};


/// @}

}// end of namespace

////////////////////////////////
//	include implementation
#include "refinement_callbacks_impl.hpp"

#endif
