// created by Sebastian Reiter
// s.b.reiter@googlemail.com
//	y10 m07 d08

#ifndef __H__LIB_GRID__REFINEMENT_CALLBACKS__
#define __H__LIB_GRID__REFINEMENT_CALLBACKS__

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
		virtual void new_vertex(VertexBase* vrt, VertexBase* parent) = 0;
		virtual void new_vertex(VertexBase* vrt, EdgeBase* parent) = 0;
		virtual void new_vertex(VertexBase* vrt, Face* parent) = 0;
		virtual void new_vertex(VertexBase* vrt, Volume* parent) = 0;
};



////////////////////////////////////////////////////////////////////////
///	calculates new positions of vertices through linear interpolation.
/**	Make sure to initialise the callback correctly. Use the same grid
 *	on which the refinement-operations will be performed. Make sure
 *	that aPos (given in the constructor) is attached to the vertices
 *	of the grid.
 *
 *	An uninitialized refinement-callback may not be used during refinement.
 */
class RefinementCallbackLinear : public IRefinementCallback
{
	public:
		RefinementCallbackLinear();
		
	///	make sure that aPos is attached to the vertices of the grid.
		RefinementCallbackLinear(Grid& grid, APosition& aPos = aPosition);
	
		virtual ~RefinementCallbackLinear();
		
		virtual void new_vertex(VertexBase* vrt, VertexBase* parent);
		virtual void new_vertex(VertexBase* vrt, EdgeBase* parent);
		virtual void new_vertex(VertexBase* vrt, Face* parent);
		virtual void new_vertex(VertexBase* vrt, Volume* parent);
		
	protected:
		Grid* 										m_pGrid;
		Grid::VertexAttachmentAccessor<APosition>	m_aaPos;
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
class RefinementCallbackEdgePlaneCut : public RefinementCallbackLinear
{
	public:
		using RefinementCallbackLinear::new_vertex;
		
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

/// @}

}// end of namespace

#endif
