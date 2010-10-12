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
		
	protected:
		typedef typename TAPosition::ValueType		pos_type;
		
		Grid* 										m_pGrid;
		Grid::VertexAttachmentAccessor<TAPosition>	m_aaPos;
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
		
	///	make sure that aPos is attached to the vertices of the grid.
		RefinementCallbackSubdivBoundary(MultiGrid& mg,
										TAPosition& aPos);
	
		virtual ~RefinementCallbackSubdivBoundary();
		
		virtual void new_vertex(VertexBase* vrt, VertexBase* parent);
		virtual void new_vertex(VertexBase* vrt, EdgeBase* parent);
		
	protected:
		virtual bool is_crease_vertex(VertexBase* vrt);
		virtual bool is_crease_edge(EdgeBase* edge);

	protected:
		MultiGrid* m_pMG;
};


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
		
	///	make sure that aPos is attached to the vertices of the grid.
		RefinementCallbackSubdivisionLoop(MultiGrid& mg,
										  TAPosition& aPos);
	
		virtual ~RefinementCallbackSubdivisionLoop();
		
		virtual void new_vertex(VertexBase* vrt, VertexBase* parent);
		virtual void new_vertex(VertexBase* vrt, EdgeBase* parent);
		virtual void new_vertex(VertexBase* vrt, Face* parent);
		
	protected:
		virtual bool is_crease_vertex(VertexBase* vrt);
		virtual bool is_crease_edge(EdgeBase* edge);

	protected:
		MultiGrid* m_pMG;
};


/// @}

}// end of namespace

////////////////////////////////
//	include implementation
#include "refinement_callbacks_impl.hpp"

#endif
